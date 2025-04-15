#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
@File    :		_encoder.py
@Time    :   	2024/06/05 16:31:25
@Author  :   	Yong Bai 
@Contact :   	baiyong at genomics.cn
@License :   	(C)Copyright 2023-2024, Yong Bai

                Licensed under the Apache License, Version 2.0 (the "License");
                you may not use this file except in compliance with the License.
                You may obtain a copy of the License at

                    http://www.apache.org/licenses/LICENSE-2.0

                Unless required by applicable law or agreed to in writing, software
                distributed under the License is distributed on an "AS IS" BASIS,
                WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
                See the License for the specific language governing permissions and
                limitations under the License.

@Desc    :   	None

"""

import math
import torch
import torch.nn as nn
import torch.nn.functional as F

try:
    from flash_attn.modules.mha import MHA
    flash_attn_available = True
except ImportError:
    import warnings
    warnings.warn("flash_attn is not available, use regular attention calculation.")
    flash_attn_available = False


class EncoderBlock(nn.Module):

    def __init__(
        self, 
        input_dim, 
        num_heads, 
        ffn_dim,
        cross_attn=False,
        dwconv=False,
        ffn_type='mlp', # MLP or SwitchMoE or  
        dropout=0.0,
        num_experts=8,
        moe_topk=4,
        use_flash_attn=False,
        device=None,
        dtype=None,
    ):
        """
        Inputs:
            input_dim - Dimensionality of the input
            num_heads - Number of heads to use in the attention block
            dim_feedforward - Dimensionality of the hidden layer in the MLP
            dropout - Dropout probability to use in the dropout layers
        """
        super().__init__()

        factory_kwargs = {'device': device, 'dtype': dtype}

        self.cross_attn = cross_attn

        # Attention layer
        self.self_attn = MHA(
            embed_dim=input_dim,
            num_heads=num_heads,
            cross_attn=cross_attn,
            dropout=dropout,
            dwconv=dwconv,
            use_flash_attn=use_flash_attn,  # if True, assert qkv.dtype in [torch.float16, torch.bfloat16]
            **factory_kwargs,
        ) if flash_attn_available else nn.MultiheadAttention(
            embed_dim=input_dim, 
            num_heads=num_heads, 
            dropout=dropout, 
            batch_first=True,
            **factory_kwargs,
        )

        if ffn_type == 'mlp':
            self.ffn = MLP(
                input_dim, 
                hidden_dim=ffn_dim, 
                output_dim=None,
                dropout=dropout, 
                multiple_of=128,
                **factory_kwargs)
        elif ffn_type == 'moe':
            self.ffn = SwitchMoE(
                input_dim, 
                ffn_dim, 
                output_dim=None,
                num_experts=num_experts, 
                k=moe_topk,
                dropout=dropout,
                **factory_kwargs
            )
        elif ffn_type == 'gated_mlp':
            self.ffn = GatedMLP(
                in_features=input_dim, 
                hidden_features=ffn_dim, 
                out_features=None, 
                return_residual=False,
                multiple_of=128, 
                **factory_kwargs
            )
        else:
            raise NotImplementedError("ffn_type must be ['mlp' | 'moe']")

        # Layers to apply in between the main layers
        self.norm1 = nn.LayerNorm([input_dim], **factory_kwargs)
        self.norm2 = nn.LayerNorm([input_dim], **factory_kwargs)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x, x_kv=None, key_padding_mask=None):
        # Attention part

        z = self.norm1(x)

        if flash_attn_available and key_padding_mask is None:
            attn_out = self.self_attn(
                z, x_kv=x_kv
            ) 
        else:
            attn_out, _ = self.self_attn(
                z, 
                z if x_kv is None else x_kv, 
                z if x_kv is None else x_kv, 
                key_padding_mask=key_padding_mask, need_weights=False
            )
        x = x + self.dropout(attn_out)

        # FFN part
        z = self.norm2(x)
        ffn_out = self.ffn(z)
        x = x + self.dropout(ffn_out)
        
        return x


class TransformerEncoder(nn.Module):

    def __init__(
        self,
        input_dim, 
        num_heads, 
        num_layers,
        ffn_dim,
        cross_attn=False,
        dwconv=False,
        ffn_type='mlp', # MLP or SwitchMoE 
        dropout=0.1,
        num_experts=8,
        moe_topk=4, 
        use_flash_attn=False,
        device=None, 
        dtype=None
    ):

        super().__init__()

        factory_kwargs = {'device': device, 'dtype': dtype}

        self.layers = nn.ModuleList(
            [
                EncoderBlock(
                    input_dim=input_dim, 
                    num_heads=num_heads,
                    ffn_dim=ffn_dim,
                    cross_attn=cross_attn,
                    ffn_type=ffn_type,
                    dwconv=dwconv,
                    dropout=dropout,
                    num_experts=num_experts,
                    moe_topk=moe_topk,
                    use_flash_attn=use_flash_attn,
                    **factory_kwargs,
                ) for _ in range(num_layers)
            ]
        )

        self.norm = nn.LayerNorm([input_dim], **factory_kwargs)
    
    def forward(self, x, x_kv=None, key_padding_mask=None):
        for l in self.layers:
            x = l(x, x_kv=x_kv, key_padding_mask=key_padding_mask)
        x = self.norm(x)
        return x

    def get_attention_maps(self, x, x_kv=None, key_padding_mask=None):
        assert not flash_attn_available, "You are using flash_attn, which does not support returning attention_weight." 
        attention_maps = []
        for l in self.layers:
            _, attn_map = l.self_attn(x, x_kv, x_kv, key_padding_mask=key_padding_mask, need_weights=True)
            attention_maps.append(attn_map)
            x = l(x)
        return attention_maps

class MLP(nn.Module):
    def __init__(
        self, 
        input_dim,
        hidden_dim=None,
        output_dim=None,
        dropout=0.1,
        multiple_of=128,
        bias1=True,
        bias2=True,
        device=None,
        dtype=None,
    ) -> None:
        super().__init__()
        factory_kwargs = {'device': device, 'dtype': dtype}
        output_dim = output_dim if output_dim is not None else input_dim
        hidden_dim = (
            hidden_dim if hidden_dim is not None else int(8 * input_dim / 3)
        )

        hidden_dim = (hidden_dim + multiple_of - 1) // multiple_of * multiple_of

        self.fc1 = nn.Linear(input_dim, hidden_dim, bias=bias1, **factory_kwargs)
        self.fc2 = nn.Linear(hidden_dim, output_dim, bias=bias2, **factory_kwargs)
        self.act = nn.GELU()
        self.dropout = nn.Dropout(dropout)

    def forward(self, x):

        x = self.act(self.fc1(x))
        x = self.dropout(x)
        return self.fc2(x)
    

# Adapted from https://github.com/Dao-AILab/flash-attention/blob/main/flash_attn/modules/mlp.py#L99C1-L136C57
class GatedMLP(nn.Module):
    def __init__(
        self,
        in_features,
        hidden_features=None,
        out_features=None,
        bias1=True,
        bias2=True,
        multiple_of=128,
        return_residual=False,
        device=None,
        dtype=None,
    ):
        factory_kwargs = {"device": device, "dtype": dtype}
        super().__init__()
        out_features = out_features if out_features is not None else in_features
        hidden_features = (
            hidden_features if hidden_features is not None else int(8 * in_features / 3)
        )
        hidden_features = (hidden_features + multiple_of - 1) // multiple_of * multiple_of
        self.return_residual = return_residual
        self.fc1 = nn.Linear(in_features, 2 * hidden_features, bias=bias1, **factory_kwargs)
        self.fc2 = nn.Linear(hidden_features, out_features, bias=bias2, **factory_kwargs)

    def forward(self, x):
        y = self.fc1(x)
        y = F.glu(y, dim=-1)
        y = self.fc2(y)
        return y if not self.return_residual else (y, x)


# Switch Transformer Mixture of Experts (MoE) Block
# NOTE: that this implementation is not optimized for speed.
# GPU memory usage keep increasing with the number of experts.
class SwitchMoE(nn.Module):
    def __init__(
        self, 
        input_dim, 
        ffn_dim, 
        output_dim=None,
        num_experts=8, 
        k=2, 
        dropout=0.1,
        device=None,
        dtype=None,
    ):
        super(SwitchMoE, self).__init__()
        factory_kwargs={'device': device, 'dtype': dtype}
        
        output_dim = output_dim if output_dim is not None else input_dim

        self.num_experts = num_experts
        self.k = k  # Number of experts to route to

        # Expert networks
        self.experts = nn.ModuleList([
            MLP(
                input_dim=input_dim, 
                hidden_dim=ffn_dim, 
                output_dim=output_dim, 
                dropout=0.0, 
                **factory_kwargs
            ) for _ in range(num_experts)
        ])

        # Gating network
        self.gate = nn.Linear(input_dim, num_experts, bias=False, **factory_kwargs)

    def forward(self, x):
        gate_scores = self.gate(x)
        top_k_gate_scores, top_k_indices = torch.topk(gate_scores, self.k, dim=-1)

        # Normalize top-k gate scores
        top_k_gate_probs = torch.softmax(top_k_gate_scores, dim=-1)

        # Initialize output
        output = torch.zeros_like(x)

        for i in range(self.k):
            selected_expert_index = top_k_indices[:, :, i]
            batch_size, seq_len = x.size(0), x.size(1)
            
            selected_experts = torch.stack([
                self.experts[selected_expert_index[b, seqlen_idx]](x[b, seqlen_idx, :]) 
                for b in range(batch_size) 
                for seqlen_idx in range(seq_len)], dim=0)
            
            selected_experts = selected_experts.view(batch_size, seq_len, -1).contiguous()
               
            output += top_k_gate_probs[:, :, i].unsqueeze(-1) * selected_experts

        return output


class OutputLayer(nn.Module):
    def __init__(
        self,
        input_dim,
        intermediate_dim=None,
        n_out=1,
        dropout=0.1,
        device=None,
        dtype=None,
    ) -> None:
        super().__init__()
        factory_kwargs = {'device': device, 'dtype': dtype}
        if intermediate_dim is None:
            intermediate_dim = input_dim

        self.out = nn.Sequential(
            # nn.LayerNorm(input_dim, **factory_kwargs), # LayerNorm is not necessary, and may cause the model to diverge (grad norm = inf)
            nn.Linear(input_dim, intermediate_dim, **factory_kwargs),
            nn.GELU(),
            # nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(intermediate_dim, n_out, **factory_kwargs),
            # nn.ReLU(),
        )
    
    def forward(self, x):
        return self.out(x)
    
    
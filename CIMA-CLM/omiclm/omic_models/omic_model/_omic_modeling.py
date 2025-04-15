#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
@File    :		_omic_modeling.py
@Time    :   	2024/06/05 11:51:16
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
import os
from typing import Dict, Optional, Union
import math
import torch
import torch.nn as nn
import torch.nn.functional as F

import time

from einops import rearrange

from transformers import PreTrainedModel, PretrainedConfig

from ._pooling import OmicInputProjection, EmbeddingPooling
from ._encoder import (
    TransformerEncoder,
    OutputLayer,
)

# to check for NaNs
# torch.autograd.set_detect_anomaly(True)

class OmicFormerConfig(PretrainedConfig):

    model_type = "omicformer"

    def __init__(
        self,
        seq_input_pooling_size: int = 501,
        
        pre_layer_type: str = 'gated_mlp', # 'gated_mlp', 'mlp
        
        ffn_type: str = 'gated_mlp', # 'moe', 'mlp', 'gated_mlp'
        hidden_dim: int = 512,
        intermediate_hidden_dim: int = 768,
        
        n_layers_encoder: int = 4,
        
        n_heads: int = 8, 
        num_experts: int = 4,
        moe_topk:int = 2,
        dropout: float = 0.1,

        initializer_range: float = 0.02,

        fusion_type: str = 'cross_attn', # 'cross_attn', 'clip_style'
        n_layers_fusion: int = 8,

        n_residuals_per_layer: int = 2,     # Change to 2 if we have MLP, otherwise 1

        out_pooling_size: int = 4,
        out_pooling_mode: str = 'adaptive',   # fisrt, last, mean, weight, adaptive, max

        n_outputs: int = 1,

        **kwargs
    ):

        self.seq_input_pooling_size = seq_input_pooling_size
        
        self.pre_layer_type = pre_layer_type

        self.ffn_type = ffn_type
        self.hidden_dim = hidden_dim
        self.intermediate_hidden_dim = intermediate_hidden_dim

        self.n_layers_encoder = n_layers_encoder

        self.n_heads = n_heads
        self.num_experts = num_experts
        self.moe_topk = moe_topk
        self.dropout = dropout

        self.initializer_range = initializer_range

        self.fusion_type = fusion_type
        self.n_layers_fusion = n_layers_fusion

        self.n_residuals_per_layer = n_residuals_per_layer

        self.out_pooling_size = out_pooling_size
        self.out_pooling_mode = out_pooling_mode

        self.n_outputs = n_outputs

        # super().__init__(architectures = ["OmicFormer"], **kwargs)
        super().__init__(**kwargs)


class OmicFormerPreTrainedModel(PreTrainedModel):
    """OmicFormerPreTrainedModel
    """
    config_class = OmicFormerConfig
    base_model_prefix = "omicformer"
    # main_input_name = "input_ids"

    def __init__(
        self, 
        config: OmicFormerConfig,
        
        seq_emb_dim=256,
        seq_emb_model = None,
        seq_emb_extraction_kwargs=None,

        scrna_emb_dim=512, 
        scrna_emb_model = None,   # for future extension
        scrna_emb_extraction_kwargs=None, # for future extenstion

        device = None,
        dtype = None,
        **kwargs
    ):
        factory_kwargs = {'device': device, 'dtype': dtype}
        super().__init__(config, **kwargs)

        self.seq_emb_model = seq_emb_model
        self.seq_emb_extraction_kwargs = seq_emb_extraction_kwargs
        self.omicformer_config = config.to_dict()
        self.seq_input_pooling_size = self.omicformer_config.pop('seq_input_pooling_size', 320)

        self.model = OmicFormer(
            seq_emb_dim=seq_emb_dim,
            scrna_emb_dim=scrna_emb_dim, 

            **self.omicformer_config,
            **factory_kwargs,
            **kwargs,
        )

    def forward(
        self,  
        cell_embedding=None,
        seq_input_ids=None, 
        seq_len=None,
        **kwargs):

        # t_start = time.time()
        seq_embedding = self._extract_embedding(
            seq_input_ids, seq_len, 
            **(self.seq_emb_extraction_kwargs if self.seq_emb_extraction_kwargs is not None else {}))
        # torch.cuda.current_stream().synchronize()  
        # t_seq_emd_end = time.time()
        # print(f"seq embedding time: {t_seq_emd_end - t_start}")
        
        out = self.model(
            seq_embedding=seq_embedding,
            cell_embedding=cell_embedding, 
            **kwargs)
        
        # return {"logits": out.squeeze(-1) if self.config.n_outputs == 1 else out}
        return {"logits": out}
    
    @torch.no_grad()
    def _extract_embedding(
        self, 
        input_ids, 
        seq_len=None,
        target='seq', # for future extension
        **kwargs
    ):
        _seq_embedding = self.seq_emb_model.extract_sample_embedding(input_ids, **kwargs)
        
        # paddeding and pooling (left padding)
        seq_embedding = torch.stack([
            (F.adaptive_avg_pool1d(
                _seq_embedding[idx, -seq_len[idx]:, :].transpose(1, 0), 
                self.seq_input_pooling_size)
            ).transpose(1, 0) 
            for idx in range(_seq_embedding.size(0))
        ])
        return seq_embedding
    
    def save_pretrained(self, save_directory,
                        congif_file_name="config.json", 
                        model_file_name="pytorch_model.bin"):
        """
        Minimal implementation of save_pretrained for MambaLMHeadModel.
        Save the model and its configuration file to a directory.
        """
        # Ensure save_directory exists
        if not os.path.exists(save_directory):
            # multiple GPUs will raise FileExistsError if exist_ok=False(default value)
            os.makedirs(save_directory, exist_ok=True)

        # Save the model's state_dict
        model_path = os.path.join(save_directory, model_file_name)
        torch.save(self.state_dict(), model_path)

        # Save the configuration of the model
        config_path = os.path.join(save_directory, congif_file_name)
        self.config.to_json_file(config_path)


class OmicFormer(nn.Module):
    """OmicFormer model
    """
    def __init__(
        self, 
        seq_emb_dim: int = 256,
        scrna_emb_dim: int = 512,
        
        pre_layer_type: str = 'gated_mlp', # 'gated_mlp', 'mlp
        
        ffn_type: str = 'gated_mlp', # 'moe', 'mlp
        hidden_dim: int = 512,
        intermediate_hidden_dim: int = 768,
        
        n_layers_encoder: int = 4,
        
        n_heads: int = 8, 
        num_experts: int = 4,
        moe_topk:int = 2,
        dropout: float = 0.1,

        initializer_range: float = 0.02,

        fusion_type: str = 'cross_attn', # 'cross_attn', 'clip_style'
        n_layers_fusion: int = 8,

        n_residuals_per_layer: int = 2,     # Change to 2 if we have MLP, otherwise 1

        out_pooling_size: int = 4,
        out_pooling_mode: str = 'adaptive',   # fisrt, last, mean, weight, adaptive, max
        n_outputs: int = 1,

        device=None,
        dtype=None,
        
        **kwargs,
        
    ) -> None:
        
        factory_kwargs = {'device': device, 'dtype': dtype}
        super().__init__()

        self.initializer_range = initializer_range
        self.n_layers_encoder = n_layers_encoder
        self.n_layers_fusion = n_layers_fusion
        self.n_residuals_per_layer = n_residuals_per_layer
        
        # seq embedding, get rid of padding token and pooling
        self.seq_input_proj = OmicInputProjection(
            emb_dim_input=seq_emb_dim,
            emb_proj_type=pre_layer_type,
            emb_proj_hidden_dim = intermediate_hidden_dim,
            hidden_dim=hidden_dim,
            **factory_kwargs
        )

        # # scrna embedding, get rid of padding token and pooling
        # # NOTE: here, we have scrna embedded, but seq embedding is too large to embed in advance.
        # # And so here we only do: removal of padding and pooling
        self.scrna_input_proj = OmicInputProjection(
            emb_dim_input=scrna_emb_dim,
            emb_proj_type=pre_layer_type,
            emb_proj_hidden_dim = intermediate_hidden_dim,
            hidden_dim=hidden_dim,
            **factory_kwargs
        )

        self.seq_emb_norm = nn.LayerNorm([hidden_dim], **factory_kwargs)
        self.scrna_emb_norm = nn.LayerNorm([hidden_dim], **factory_kwargs)

        self.seq_encoder = TransformerEncoder(
            hidden_dim, 
            n_heads,
            n_layers_encoder,
            ffn_dim = intermediate_hidden_dim,
            cross_attn=False,
            dwconv=True,
            ffn_type=ffn_type,
            dropout=dropout,
            num_experts=num_experts,
            moe_topk=moe_topk,
            **factory_kwargs,
        )

        self.scrna_encoder = TransformerEncoder(
            hidden_dim, 
            n_heads,
            n_layers_encoder,
            ffn_dim = intermediate_hidden_dim,
            cross_attn=False,
            dwconv=True,
            ffn_type=ffn_type,
            dropout=dropout,
            num_experts=num_experts,
            moe_topk=moe_topk,
            **factory_kwargs,
        )

        # we fuse the seq and scrna using cross attention
        if fusion_type == 'cross_attn':
            self.fusion = TransformerEncoder(
                hidden_dim, 
                n_heads,
                n_layers_fusion,
                ffn_dim = intermediate_hidden_dim,
                cross_attn=True,
                dwconv=False,
                ffn_type=ffn_type,
                dropout=dropout,
                num_experts=num_experts,
                moe_topk=moe_topk,
                **factory_kwargs,
            )
        elif fusion_type == 'clip_style':
            raise NotImplementedError("clip_style fusion is not implemented yet.")

        self.out_pooling = EmbeddingPooling(
            pooling_seq_len=out_pooling_size,
            pooling_mode=out_pooling_mode,  
            **factory_kwargs
        )
        
        # prediction layer
        self.out_proj = OutputLayer(
            input_dim=hidden_dim * out_pooling_size,
            intermediate_dim=intermediate_hidden_dim,
            n_out=n_outputs,
            dropout=dropout,
            **factory_kwargs)
        
        self.out_pooling_size = out_pooling_size

        
        self.apply(self._init_weights)

    # https://github.com/huggingface/transformers/blob/c28d04e9e252a1a099944e325685f14d242ecdcd/src/transformers/models/gpt2/modeling_gpt2.py#L454
    def _init_weights(self, module):
        """Initialize the weights."""
        if isinstance(module, nn.Linear):
            # Slightly different from the TF version which uses truncated_normal for initialization
            # cf https://github.com/pytorch/pytorch/pull/5617
            module.weight.data.normal_(mean=0.0, std=self.initializer_range)
            if module.bias is not None:
                if not getattr(module.bias, "_no_reinit", False):
                    nn.init.zeros_(module.bias)
        elif isinstance(module, nn.LayerNorm):
            module.bias.data.zero_()
            module.weight.data.fill_(1.0)

        # # Reinitialize selected weights subject to the OpenAI GPT-2 Paper Scheme:
        # #   > A modified initialization which accounts for the accumulation on the residual path with model depth. Scale
        # #   > the weights of residual layers at initialization by a factor of 1/√N where N is the # of residual layers.
        # #   >   -- GPT-2 :: https://openai.com/blog/better-language-models/
        # #
        # # Reference (Megatron-LM): https://github.com/NVIDIA/Megatron-LM/blob/main/megatron/model/gpt_model.py
        # for name, p in module.named_parameters():
        #     if name == "c_proj.weight":
        #         # Special Scaled Initialization --> There are 2 Layer Norms per Transformer Block
        #         p.data.normal_(
        #             mean=0.0, 
        #             std=(self.initializer_range / math.sqrt(
        #                 self.n_residuals_per_layer * (self.n_layers_encoder + self.n_layers_fusion)
        #                 ))
        #         )
        #     if name in ["out_proj.weight", "fc2.weight"]:
        #         # Special Scaled Initialization --> There are 2 Layer Norms per Transformer Block
        #         # Following Pytorch init, except scale by 1/sqrt(2 * n_layer)
        #         # We need to reinit p since this code could be called multiple times
        #         # Having just p *= scale would repeatedly scale it down
        #         nn.init.kaiming_uniform_(p, a=math.sqrt(5))
        #         with torch.no_grad():
        #             p /= math.sqrt(
        #                 self.n_residuals_per_layer * (self.n_layers_encoder + self.n_layers_fusion)
        #             )

    def forward(
        self,  
        seq_embedding=None, 
        cell_embedding=None, 
        # seq_input_ids=None,
    ):
        # seq_embedding: (batch_size, seq_len, seq_emb_dim)
        # cell_embedding: (batch_size, num_cells, scrna_emb_dim)
        # NOTE:
        # - we have embedded the scRNA cells. However, due to large amount of seq data,
        #   we embed seq data here on the fly.
        # - `seq_len` may have different length across individuals and site, 
        # - we have padded them the same length. See `OmicDataCollator` So does `num_cells`
        # 
        # dict_keys(['seq_input_ids', 'peak_value', 'cell_embedding', 'seq_len'])

        # t_start = time.time()
        seq_emb = self.seq_input_proj(seq_embedding) # seq_emb: (batch_size, seq_input_pooling_size, seq_emb_dim)
        # torch.cuda.current_stream().synchronize()  
        # t_seq_emd_end = time.time()
        # print(f"seq embedding time: {t_seq_emd_end - t_start}")

        scrna_emb = self.scrna_input_proj(cell_embedding) # scrna_emb: (batch_size, scrna_input_pooling_size, scrna_emb_dim)

        # # seq_emb encoder
        seq_emb = self.seq_emb_norm(seq_emb)
        seq_emb = self.seq_encoder(seq_emb)  # seq_emb: (batch_size, seq_input_pooling_size, hidden_dim)

        # # scrna_emb encoder
        scrna_emb = self.scrna_emb_norm(scrna_emb)
        scrna_emb = self.scrna_encoder(scrna_emb)  # scrna_emb: (batch_size, scrna_input_pooling_size, hidden_dim)

        # # cross attention
        fusion_emb = self.fusion(seq_emb, scrna_emb)

        fusion_emb = self.out_pooling(fusion_emb) # fusion_emb: (batch_size, out_pooling_size, hidden_dim)
        
        out = rearrange(fusion_emb, 'b n d -> b (n d)', n=fusion_emb.size(1))  # concat
        out = self.out_proj(out)  # out: (batch_size, n_outputs)

        return out
    

class OmicClipFormer(nn.Module):
    """OmicClipFormer model
    """
    def __init__(
        self, 
        seq_model_name: str = None,
        seq_pretrained_model_name_or_path: str = None,
        seq_emb_dim: int = None,
        seq_emb_extraction_kwargs: Dict = None,

        scrna_model_name: str = None,
        scrna_pretrained_model_name_or_path: str = None,
        scrna_emb_dim: int = None,
        scrna_emb_extraction_kwargs: Dict = None,

        seq_emb_input: bool = False,
        scrna_emb_input: bool = True,

        seq_input_pooling_size: int = 100,
        scrna_input_pooling_size: int = 100,
        seq_input_pooling_mode: str = 'adaptive',   # fisrt, last, mean, weight, adaptive, max
        scrna_input_pooling_mode: str = 'adaptive',
        
        pre_layer_type: str = 'gated_mlp', # 'gated_mlp', 'mlp
        
        ffn_type: str = 'moe', # 'moe', 'mlp
        hidden_dim: int = 768,
        intermediate_hidden_dim: int = 1024,
        
        n_layers_encoder: int = 12,
        
        n_heads: int = 8, 
        num_experts: int = 8,
        moe_topk:int = 4,
        dropout: float = 0.1,

        initializer_range: float = 0.02,

        fusion_type: str = 'cross_attn', # 'cross_attn', 'clip_style'
        n_layers_fusion: int = 12,

        n_residuals_per_layer: int = 2,     # Change to 2 if we have MLP, otherwise 1

        out_pooling_size: int = 4,
        out_pooling_mode: str = 'adaptive',   # fisrt, last, mean, weight, adaptive, max
        n_outputs: int = 1,

        device=None,
        dtype=None,
        
        **kwargs,
        
    ) -> None:
        
        factory_kwargs = {'device': device, 'dtype': dtype}
        super().__init__()

        self.initializer_range = initializer_range
        self.n_layers_encoder = n_layers_encoder
        self.n_layers_fusion = n_layers_fusion
        self.n_residuals_per_layer = n_residuals_per_layer
        
        # seq embedding, get rid of padding token and pooling
        self.seq_embedding = OmicEmbedding(
            model_name=seq_model_name,
            pretrained_model_name_or_path=seq_pretrained_model_name_or_path,
            emb_dim=seq_emb_dim,
            emb_straction_kwargs=seq_emb_extraction_kwargs,
            emb_input=seq_emb_input,
            emb_pooling_len=seq_input_pooling_size,
            emb_pooling_mode=seq_input_pooling_mode,
            emb_proj_type=pre_layer_type,
            emb_proj_hidden_dim = intermediate_hidden_dim,
            hidden_dim=hidden_dim,
            **factory_kwargs
        )

        # scrna embedding, get rid of padding token and pooling
        # NOTE: here, we have scrna embedded, but seq embedding is too large to embed in advance.
        # And so here we only do: removal of padding and pooling
        self.scrna_embedding = OmicEmbedding(
            model_name=scrna_model_name,
            pretrained_model_name_or_path=scrna_pretrained_model_name_or_path,
            emb_dim=scrna_emb_dim,
            emb_straction_kwargs=scrna_emb_extraction_kwargs,
            emb_input=scrna_emb_input,
            emb_pooling_len=scrna_input_pooling_size,
            emb_pooling_mode=scrna_input_pooling_mode,
            emb_proj_type=pre_layer_type,
            emb_proj_hidden_dim = intermediate_hidden_dim,
            hidden_dim=hidden_dim,
            **factory_kwargs
        )
        self.seq_emb_norm = nn.LayerNorm([hidden_dim], **factory_kwargs)
        self.scrna_emb_norm = nn.LayerNorm([hidden_dim], **factory_kwargs)

        self.seq_encoder = TransformerEncoder(
            hidden_dim, 
            n_heads,
            n_layers_encoder,
            ffn_dim = intermediate_hidden_dim,
            cross_attn=False,
            dwconv=True,
            ffn_type=ffn_type,
            dropout=dropout,
            num_experts=num_experts,
            moe_topk=moe_topk,
            **factory_kwargs,
        )

        self.scrna_encoder = TransformerEncoder(
            hidden_dim, 
            n_heads,
            n_layers_encoder,
            ffn_dim = intermediate_hidden_dim,
            cross_attn=False,
            dwconv=True,
            ffn_type=ffn_type,
            dropout=dropout,
            num_experts=num_experts,
            moe_topk=moe_topk,
            **factory_kwargs,
        )

        # we fuse the seq and scrna using cross attention
        if fusion_type == 'cross_attn':
            self.fusion = TransformerEncoder(
                hidden_dim, 
                n_heads,
                n_layers_fusion,
                ffn_dim = intermediate_hidden_dim,
                cross_attn=True,
                dwconv=False,
                ffn_type=ffn_type,
                dropout=dropout,
                num_experts=num_experts,
                moe_topk=moe_topk,
                **factory_kwargs,
            )
        elif fusion_type == 'clip_style':
            raise NotImplementedError("clip_style fusion is not implemented yet.")

        self.out_pooling = EmbeddingPooling(
            d_model=hidden_dim,
            pooling_seq_len=out_pooling_size,
            pooling_mode=out_pooling_mode,  
            **factory_kwargs
        )
        
        # prediction layer
        self.out_proj = OutputLayer(
            input_dim=hidden_dim * out_pooling_size,
            intermediate_dim=intermediate_hidden_dim,
            n_out=n_outputs,
            dropout=dropout,
            **factory_kwargs)
        self.out_pooling_size = out_pooling_size

        
        self.apply(self._init_weights)

    # https://github.com/huggingface/transformers/blob/c28d04e9e252a1a099944e325685f14d242ecdcd/src/transformers/models/gpt2/modeling_gpt2.py#L454
    def _init_weights(self, module):
        """Initialize the weights."""
        if isinstance(module, nn.Linear):
            # Slightly different from the TF version which uses truncated_normal for initialization
            # cf https://github.com/pytorch/pytorch/pull/5617
            module.weight.data.normal_(mean=0.0, std=self.initializer_range)
            if module.bias is not None:
                if not getattr(module.bias, "_no_reinit", False):
                    nn.init.zeros_(module.bias)
        elif isinstance(module, nn.LayerNorm):
            module.bias.data.zero_()
            module.weight.data.fill_(1.0)

        # # Reinitialize selected weights subject to the OpenAI GPT-2 Paper Scheme:
        # #   > A modified initialization which accounts for the accumulation on the residual path with model depth. Scale
        # #   > the weights of residual layers at initialization by a factor of 1/√N where N is the # of residual layers.
        # #   >   -- GPT-2 :: https://openai.com/blog/better-language-models/
        # #
        # # Reference (Megatron-LM): https://github.com/NVIDIA/Megatron-LM/blob/main/megatron/model/gpt_model.py
        # for name, p in module.named_parameters():
        #     if name == "c_proj.weight":
        #         # Special Scaled Initialization --> There are 2 Layer Norms per Transformer Block
        #         p.data.normal_(
        #             mean=0.0, 
        #             std=(self.initializer_range / math.sqrt(
        #                 self.n_residuals_per_layer * (self.n_layers_encoder + self.n_layers_fusion)
        #                 ))
        #         )
        #     if name in ["out_proj.weight", "fc2.weight"]:
        #         # Special Scaled Initialization --> There are 2 Layer Norms per Transformer Block
        #         # Following Pytorch init, except scale by 1/sqrt(2 * n_layer)
        #         # We need to reinit p since this code could be called multiple times
        #         # Having just p *= scale would repeatedly scale it down
        #         nn.init.kaiming_uniform_(p, a=math.sqrt(5))
        #         with torch.no_grad():
        #             p /= math.sqrt(
        #                 self.n_residuals_per_layer * (self.n_layers_encoder + self.n_layers_fusion)
        #             )

    def forward(
        self, 
        seq_input_ids=None, 
        cell_input_ids=None, 
        seq_embedding=None, 
        cell_embedding=None, 
        seq_len=None, 
    ):
        # seq_embedding: (batch_size, seq_len, seq_emb_dim)
        # cell_embedding: (batch_size, num_cells, scrna_emb_dim)
        # NOTE:
        # - we have embedded the scRNA cells. However, due to large amount of seq data,
        #   we embed seq data here on the fly.
        # - `seq_len` may have different length across individuals and site, 
        # - we have padded them the same length. See `OmicDataCollator` So does `num_cells`
        # 
        # dict_keys(['seq_input_ids', 'peak_value', 'cell_embedding', 'seq_len'])
    
        seq_emb = self.seq_embedding(
            input_ids=seq_input_ids,
            input_embedding=seq_embedding,
            seq_len=seq_len
        ) # seq_emb: (batch_size, seq_emb_dim)

        # print(f"seq embedded: {seq_emb.size()}")

        scrna_emb = self.scrna_embedding(
            input_ids=cell_input_ids,
            input_embedding=cell_embedding,
            seq_len=None
        ) # scrna_emb: (batch_size,scrna_emb_dim)

        seq_emb = self.seq_emb_norm(seq_emb)
        # scrna_emb encoder
        scrna_emb = self.scrna_emb_norm(scrna_emb)



        # print(f"scrna encoded: {scrna_emb.size()}")

        # cross attention
        fusion_emb = self.fusion(seq_emb, scrna_emb)

        fusion_emb = self.out_pooling(fusion_emb) # fusion_emb: (batch_size, out_pooling_size, hidden_dim)

        out = rearrange(fusion_emb, 'b n d -> b (n d)', n=fusion_emb.size(1))  # concat

        # print(f"out size {out.size()}")

        return self.out_proj(out)



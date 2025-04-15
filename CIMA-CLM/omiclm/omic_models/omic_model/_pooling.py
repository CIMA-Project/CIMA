#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
@File    :		_aggr.py
@Time    :   	2024/06/04 09:20:42
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

from typing import Dict, Optional, Union
import torch
import torch.nn as nn
import torch.nn.functional as F
from einops import rearrange

import time

from ._encoder import GatedMLP, MLP

class WeightedAvgPooling(nn.Module):
    def __init__(
        self, 
        d_model,
        device=None,
        dtype=None
    ) -> None:

        factory_kwargs = {'device': device, 'dtype': dtype}
        super().__init__()

        self.weights = nn.Parameter(torch.randn(d_model, **factory_kwargs))

    def forward(self, x):
        # x: (batch_size, seq_len, d_model)
        attn_scores = torch.einsum('bij,j->bi', x, self.weights)    # attn_scores: (batch_size, seq_len)
        attn_weights = torch.softmax(attn_scores, dim=-1)
        weighted_sum = torch.einsum("bij,bi->bj", x, attn_weights)  # weighted_sum: (batch_size, d_model)

        return weighted_sum


class EmbeddingPooling(nn.Module):
    def __init__(
        self, 
        pooling_seq_len: Optional[int]=100,
        pooling_mode: str='adaptive',    #  fisrt, last, mean, weight, adaptive
        device=None,
        dtype=None,
    ):
        """ Obtain the `pooling_seq_len` embedding vectors from the input embedding
        metrix using the given `mode`.

        input embedding metrix: x, with shape (batch_size, seq_len, d_model)
            e.g, embedding metrix for single cell: (batch_size, num_cells, d_model) 

        Parameters
        ----------
        d_model : int
            dimension of embedding
        pooling_seq_len : int, optional
            the length of output sequence, by default 1
        pooling_mode : str, optional
            the strategy to grab `pooling_seq_len` embedding vectors, by default 'avg'.
                `first`: grab the first `pooling_seq_len` embedding vectors,
                `last`:  grab the last `pooling_seq_len` embedding vectors,
                `weighted`: grab the weighted sum of the entire embedding metrix,
                    NOTE: `weighted` or `avg` only return one weighted sum embedding vector. (ie., `pooling_seq_len=1`)
                `mean`: grab the cummulative mean of embedding metrix,
                `adaptive`: grab the `pooling_seq_len` embedding vectors using adaptive average pooling.

        device : _type_, optional
            _description_, by default None
        dtype : _type_, optional
            _description_, by default None
        """
        factory_kwargs = {'device': device, 'dtype': dtype}
        super().__init__()

        self.pooling_seq_len = pooling_seq_len
        self.pooling_mode = pooling_mode

        self.ada_pool = nn.AdaptiveAvgPool1d(pooling_seq_len)

        # if self.pooling_mode == 'weighted':
        #     self.weight_pool = WeightedAvgPooling(d_model=d_model, **factory_kwargs)

    def forward(self, x):

        squenze = False
        if len(x.size()) != 3 :
            x = x.unsqueeze(0)  # x: (batch_size, seq_len, d_model)
            squenze=True
        
        # adapted from HyenaDNA SequenceDecoder
        if self.pooling_mode == "first":
            retrive = lambda x: x[..., :self.pooling_seq_len, :]
        elif self.pooling_mode == "last":
            retrive = lambda x: x[..., -self.pooling_seq_len:, :]
        elif self.pooling_mode == 'mean':
            retrive = lambda x: (
                torch.cumsum(x, dim=-2)
                / torch.arange(
                    1, 1 + x.size(-2), device=x.device, dtype=x.dtype
                ).unsqueeze(-1)
            )[..., -self.pooling_seq_len:, :]
        elif self.pooling_mode == 'adaptive':
            retrive = lambda x: rearrange(
            self.ada_pool(rearrange(x, 'b l d -> b d l')), 
            'b d l-> b l d'
        )
        # elif self.pooling_mode == 'weighted':
        #     retrive = lambda x: self.weight_pool(x).unsqueeze(-2)
        elif self.pooling_mode == 'max':
            retrive = lambda x: torch.max(x, dim=-2, keepdim=True).values
        else:
            raise NotImplementedError("mode must be ['last' | 'first' | 'mean' | 'adaptive' | 'weighted']")
        
        x = retrive(x)

        if squenze:
            x = x.squeeze(0)
        return x
    

class OmicInputProjection(nn.Module):
    def __init__(
        self, 
        emb_dim_input: int = 256,    # embedding dim by the pretrained embedding model

        emb_proj_type: str = 'gated_mlp', # 'gated_mlp', 'mlp',
        emb_proj_hidden_dim: int = 512,
        emb_proj_dropout: float = 0.0,
        hidden_dim: int = 512,  # output dim of the embedding projection

        device=None,
        dtype=None, 
        **kwargs
    ):
        factory_kwargs = {'device': device, 'dtype': dtype}
        super().__init__(**kwargs)

        self.pre_proj = nn.Sequential(
            nn.LayerNorm([emb_dim_input], **factory_kwargs), 
            (MLP(
                emb_dim_input, 
                hidden_dim=emb_proj_hidden_dim, 
                output_dim=hidden_dim, 
                dropout=emb_proj_dropout,
                **factory_kwargs) 
             if emb_proj_type == 'mlp' 
             else GatedMLP(
                emb_dim_input, 
                hidden_features=emb_proj_hidden_dim,
                out_features=hidden_dim,
                **factory_kwargs)
            )
        )

    def forward(self, input_embedding):
        embeddings = self.pre_proj(input_embedding)
        return embeddings
    
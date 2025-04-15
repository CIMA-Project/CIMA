#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
@File    :		_geneformer_model.py
@Time    :   	2024/05/26 10:41:20
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
import logging
from typing import Any, List, Optional
import numpy as np
import torch
import torch.nn.functional as F
from tqdm.auto import trange
from anndata import AnnData
import scanpy as sc

from transformers import (
    BertForMaskedLM,
    PreTrainedModel,
    BertConfig,
)
from datasets import Dataset

logger = logging.getLogger(__name__)

def pad_tensor(t: torch.Tensor,
               max_size: int,
               pad_token_id: int = 0) -> torch.Tensor:
    """
    Pad a tensor to a max size
    """
    
    return F.pad(t, pad = (0, max_size - t.numel()), 
                 mode = 'constant', value = pad_token_id)

class GeneformerModelWrapper(PreTrainedModel):
    config_class = BertConfig

    def __init__(
        self,
        config=None,
        model: BertForMaskedLM = None,
        device=None,
        **kwargs,
    ) -> None:
        
        super().__init__(config, **kwargs)

        self.model = model
        self.model.to(device)
        
        self.config = self.model.config
    
    # adapted from https://github.com/microsoft/zero-shot-scfoundation/blob/main/sc_foundation_evals/geneformer_forward.py#L276
    def extract_sample_embedding(
        self,
        tokenized_dataset: Dataset,
        adata: AnnData,
        batch_size: int = 48,
        embedding_key: str = "cell_embedding",
        layer: int=-2,
        pad_token_id: int = 0,
        obs_to_save: Optional[List[str]] = 'sample',
        return_new_adata=True, 
        **kwargs,
    ):
        n_layers = self.model.config.num_hidden_layers
        if layer >= n_layers or layer < -n_layers:
            raise ValueError(f"Layer {layer} is not valid. There are only {n_layers} "
                   f"Acceptable values are between {-n_layers} (if counting "
                   f"forwards) and {n_layers - 1} (if counting backwards)")
       
        device = next(self.model.parameters()).device
        self.model.eval()
        
        data_size = len(tokenized_dataset)
        cell_embs_list = []
        rankings_list = []
        for i in trange(
            0, data_size, batch_size,
            desc = "Geneformer (extracting embeddings)"):

            max_range = min(i + batch_size, data_size)
            batch_dataset = tokenized_dataset.select(range(i, max_range))
            batch_dataset.set_format(type="torch")

            org_len = batch_dataset["length"].clone().to(device)

            batch, attn_mask = self._extend_batch(batch_dataset, pad_token_id=pad_token_id)
            
            model_output = self._pass_batch(batch, attention_mask = attn_mask)

            embs = model_output.hidden_states[layer]

            cell_embs = self.mean_nonpadding_embs(embs, org_len)

            # add cell embeddings to the list
            cell_embs_list.extend(cell_embs.detach().cpu().numpy())

            # now, get the ranking reconstruction
            out_rankings = (model_output.logits
                            .argmax(axis=-1)
                            .detach().cpu().numpy())
            
            # save the rankings with the original order
            rankings_list.extend(out_rankings)
            
            torch.cuda.empty_cache()
            del model_output
            del batch
            del attn_mask
            del embs
            del cell_embs

        cell_embeddings = np.array(cell_embs_list)

        output_rankings = rankings_list
        input_rankings = [np.array(item) 
                          for item 
                          in tokenized_dataset['input_ids']]

        adata.obsm["output_rankings"] = output_rankings
        adata.obsm["input_rankings"] = input_rankings

        # add embeddings to adata
        adata.obsm[embedding_key] = cell_embeddings

        if isinstance(obs_to_save, str):
            assert obs_to_save in adata.obs, f"obs_to_save {obs_to_save} not in adata.obs"
            obs_to_save = [obs_to_save]

        if return_new_adata:
            obs_df = adata.obs[obs_to_save] if obs_to_save is not None else None
            return sc.AnnData(X=cell_embeddings, obs=obs_df, dtype="float32")

        return adata


    def _extend_batch(self,
                      batch_dataset: Dataset,
                      return_attention_mask: bool = True,
                      pad_token_id: int = 0):
        
        max_size = max(batch_dataset['length'])
        
        batch_ = [pad_tensor(x, max_size, pad_token_id) 
                  for x in batch_dataset['input_ids']]
        
        device = next(self.model.parameters()).device

        batch_ = torch.stack(batch_).to(device)

        if return_attention_mask:
            mask_ = [[1] * l + [0] * (max_size - l) 
                     for l in batch_dataset['length']]
            mask_ = torch.tensor(mask_).to(device)
            return batch_, mask_
            
        return batch_
    
    def _pass_batch(self, 
                    batch_ids: torch.Tensor, 
                    attention_mask: torch.Tensor,
                    **kwargs) -> torch.Tensor:
        # make sure that batch and attn_mask on the same device
        device = next(self.model.parameters()).device

        batch_ids = batch_ids.to(device)
        attn_mask = attention_mask.to(device)

        # with torch.no_grad():
        with torch.inference_mode():
            outputs = self.model(input_ids = batch_ids,
                                 attention_mask = attn_mask,
                                 **kwargs)
        
        return outputs
    
    # get cell embeddings excluding padding
    def mean_nonpadding_embs(self, embs, original_lens):
        # mask based on padding lengths
        # mask = torch.arange(embs.size(1)).unsqueeze(0).to("cuda") < original_lens.unsqueeze(1)
        device = next(self.model.parameters()).device
        mask = torch.arange(embs.size(1)).unsqueeze(0).to(device) < original_lens.unsqueeze(1)
        # extend mask dimensions to match the embeddings tensor
        mask = mask.unsqueeze(2).expand_as(embs)

        # use the mask to zero out the embeddings in padded areas
        masked_embs = embs * mask.float()

        # sum and divide by the lengths to get the mean of non-padding embs
        mean_embs = masked_embs.sum(1) / original_lens.view(-1, 1).float()
        return mean_embs

    @classmethod
    def from_pretrained(
        cls,
        pretrained_model_name_or_path: str,
        output_attentions=False,
        output_hidden_states=True,
        local_files_only=True,
        device=None,
        **kwargs,
    ):
        geneforer_model = BertForMaskedLM.from_pretrained(
            pretrained_model_name_or_path,
            output_attentions=output_attentions,
            output_hidden_states=output_hidden_states,
            local_files_only=local_files_only,
            **kwargs
        )

        model = cls(BertConfig(), model=geneforer_model, device=device)
        return model
    
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

        



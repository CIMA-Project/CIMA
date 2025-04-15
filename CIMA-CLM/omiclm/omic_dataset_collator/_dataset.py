#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
@File    :		_dataset.py
@Time    :   	2024/05/20 09:45:01
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
from typing import Any, Dict, List, Mapping, NewType, Optional, Sequence, Union
from dataclasses import dataclass
import logging
from scipy.sparse import issparse
import itertools

import torch
import torch.nn.functional as F
import torch.utils.data
import datasets
import copy
import scanpy as sc

from anndata import AnnData

import time

logger = logging.getLogger(__name__)

InputDataClass = NewType("InputDataClass", Any)


class OmicIterableDataset(torch.utils.data.IterableDataset):
    """IterableDataset that combines seq sequence, ATAC peaks and RNA. These data is undergoing process in order to feed LLM.
    """
    def __init__(
        self,
        seq_tokenized_iter_dataset: datasets.IterableDataset,
        scrna_embedding_data,
    ) -> None:
        """OmicDataset derived from torch.utils.data.Dataset 

        Parameters
        ----------
        seq_tokenized_iter_dataset : datasets.IterableDataset
            The tokenized seq iterable dataset.

            >>> Example:
            IterableDataset({
                features: ['sample', 'pos', 'peak_value', 'input_ids', 'record_id'],
                n_shards: 500
            })

        """

        super().__init__()

        self.seq_token_ds = seq_tokenized_iter_dataset
        self.scrna_embedding_data = scrna_embedding_data
    
    def _process_dataset(self, i_seq_token_ds):

        sample_id = i_seq_token_ds['sample']

        if 'sample' not in self.scrna_embedding_data.obs.keys():
            raise ValueError("The sample is not found in the scRNA embedding adata.obs.keys().")

        cell_adata = self.scrna_embedding_data[
            self.scrna_embedding_data.obs['sample'] == sample_id]
        
        if cell_adata.shape[0] == 0:
            # raise ValueError(f"Not found scRNA embedding for sample {sample_id}")
            logger.warning(f"Not found scRNA embedding for sample {sample_id}")
            return None
        
        _cell_embedding = cell_adata.X.tolist()
        # _cell_ids = cell_adata.obs_names.to_list()

        # print(i_seq_token_ds["input_ids"])

        return dict(
            seq_input_ids=i_seq_token_ds["input_ids"],
            peak_value=i_seq_token_ds["peak_value"],
            cell_embedding =_cell_embedding, # 2d list
            # # the followings are not used in the model forward, just for recording the sample and id.
            # sample=sample_id,
            # pos=i_seq_token_ds["pos"] if "pos" in i_seq_token_ds else None,
            # cell_id = _cell_ids  # list
            sample_record_id=_seq_ds["record_id"],
        )
    
    def __iter__(self):

        worker_info = torch.utils.data.get_worker_info()

        worker_id = worker_info.id if worker_info is not None else 0
        worker_total_num = worker_info.num_workers if worker_info is not None else 1

        # Map each element using the self._process_data
        mapped_itr = map(self._process_dataset, self.seq_token_ds)

        # Add multiworker functionality
        mapped_itr = itertools.islice(mapped_itr, worker_id, None, worker_total_num)

        return mapped_itr

        # return iter(map(self._process_dataset, self.seq_token_ds))
    

class OmicDataset(torch.utils.data.Dataset):
    """Dataset that combines seq sequence, ATAC peaks and RNA. These data is undergoing process in order to feed LLM.
    """
    def __init__(
        self,
        seq_tokenized_dataset: datasets.Dataset,
        scrna_embedding_data: AnnData,
        peak_value_feature:str="log1p_norm_peak_value",
    ) -> None:
        """OmicDataset derived from torch.utils.data.Dataset 

        Parameters
        ----------
        seq_tokenized_dataset : datasets.Dataset
            The tokenized seq dataset.

            >>> Example:
            Dataset({
                features: ['sample', 'pos', 'seq', 'peak_value', 'input_ids', 'attention_mask', 'record_id'],
                num_rows: 111750940
            })

        scrna_embedding_data : AnnData
            The scRNA embedding data. saved in the AnnData format.
        
        peak_value_feature : str, default="log1p_norm_peak_value"
            The feature name of the peak value in the seq_tokenized_dataset.
        

        """

        super().__init__()

        self.seq_token_ds = seq_tokenized_dataset
        self.scrna_embedding_data = scrna_embedding_data
        self.peak_value_feature = peak_value_feature

    def __len__(self):
        return len(self.seq_token_ds)
    
    def __getitem__(self, index):

        _seq_ds = self.seq_token_ds[index]

        _seq_sampleid = _seq_ds["sample"]
        # pos = _seq_ds["pos"]
        _seq_peak_value = _seq_ds[self.peak_value_feature]
        cell_adata = self.scrna_embedding_data[
            self.scrna_embedding_data.obs['sample'] == _seq_sampleid]
        # 随机选择 10,000 条数据
        #cell_adata = sc.pp.subsample(self.scrna_embedding_data, n_obs=10000, copy=True)
        # cell_adata = self.scrna_embedding_data[self.scrna_embedding_data.obs['sample'] == 'hg19']
        if cell_adata.shape[0] == 0:
            # raise ValueError(f"Not found scRNA embedding for sample {_seq_sampleid}")
            # logger.warning(f"Not found scRNA embedding for sample {_seq_sampleid}")
            return None
        
        _cell_embedding = cell_adata.X # cell_adata.X.A.tolist() if issparse(cell_adata.X) else cell_adata.X.tolist() 
        # get corresponding scRNA embedding according the `_seq_sampleid`
        # _scrna_ds = self.scrna_emb_ds.filter(
        #         lambda batch: [example[self.scrna_sampleid_name] == _seq_sampleid for example in batch],
        #         batched=True,
        #         num_proc=self.filter_num_proc,
        #     )

        # for speed up. See https://github.com/huggingface/datasets/issues/5498
        # self.scrna_emb_ds = self.scrna_emb_ds.with_format("arrow")
        # _scrna_ds = self.scrna_emb_ds.filter(
        #         lambda table: table[table[self.scrna_sampleid_name] == _seq_sampleid],
        #         batched=True,
        #         num_proc=self.filter_num_proc,
        #         keep_in_memory=True,
        #         load_from_cache_file=False,
        #     )
        # self.scrna_emb_ds = self.scrna_emb_ds.with_format(None)
        # _scrna_ds = _scrna_ds.with_format(None)

        return dict(
            seq_input_ids=_seq_ds["input_ids"],
            # peak_value=_seq_ds["peak_value"],
            peak_value=_seq_peak_value,  # scalar
            cell_embedding =_cell_embedding, # 2d list
            #sample_record_id=_seq_ds["record_id"],
            # # >>> the followings are not used in the model forward, just for recording the sample and id.
            # sample=_seq_sampleid,
            # pos=_seq_ds["pos"] if "pos" in _seq_ds else None,
        )

@dataclass
class OmicDataCollator:
    """We only padded on the right side for seq data.
    
    """
    def __init__(
        self, 
        seq_pad_token_id:int=4,         
        sc_max_cnt: int=512, 
        seq_max_len: int=512,

        seq_emb_model=None,

    ): 
        self.seq_pad_token_id = seq_pad_token_id
        self.sc_max_cnt = sc_max_cnt
        self.seq_max_len = seq_max_len

        if seq_emb_model is not None:
            self.seq_emb_model = seq_emb_model

    def __call__(self, examples: List[InputDataClass]):

        # # Handle example is None
        # examples = [example for example in examples if example is not None]
        
        if not isinstance(examples[0], Mapping):
            examples = [vars(example) for example in examples]
        first = examples[0]
        batch = {}

        # get the max length in the batch
        # logger.info(f"Start to collate the batch with {len(examples)} examples")

        _seq_len = [len(example['seq_input_ids']) if len(example['seq_input_ids']) <= self.seq_max_len else self.seq_max_len for example in examples]
        #_seq_len = [len(example['seq_input_ids']) for example in examples]
        # max_seq_len = max(_seq_len)
        max_seq_len = self.seq_max_len
        min_cell_cnt = self.sc_max_cnt

        # if isinstance(self.sc_max_cnt, str) and self.sc_max_cnt == 'auto':
        #     min_cell_cnt = max([len(example['cell_embedding']) for example in examples]) 
        # elif isinstance(self.sc_max_cnt, int):
        #     min_cell_cnt = self.sc_max_cnt
        # else:
        #     raise ValueError(f"{self.sc_max_cnt} is not supported, using ['auto' | int]")

        # if self.padding_side == "right": 
        #     batch['input_ids'] = torch.nn.utils.rnn.pad_sequence(
        #         [torch.tensor(example["input_ids"]) for example in examples],
        #         padding_value=self.seq_pad_token_id,
        #         batch_first=True,
        #     )
        # else: # self.padding_side == "left":
        #     batch['input_ids'] = torch.nn.utils.rnn.pad_sequence(
        #         [torch.tensor(example["input_ids"][::-1]) for example in examples ],  # reverse the list and create tensors 
        #         padding_value=self.seq_pad_token_id,
        #         batch_first=True
        #     ).flip(dims=[-1])                    # reverse/flip the padded tensor in first dimension
        
        # Handling of all other possible keys.
        # Again, we will use the first element to figure out which key/values are not None for this model.
        for k, v in first.items():
            if isinstance(v, list) and isinstance(v[0], str):
                continue
            if v is not None and not isinstance(v, str):
                if k == "seq_input_ids":
                    # add pad_token_id on right side
                    # batch[k] = torch.tensor([ f[k] + (max_seq_len - len(f[k])) * [self.seq_pad_token_id] if len(f[k]) < max_seq_len else f[k][:max_seq_len] 
                    #                         for f in examples])
                    
                    # add pad_token_id on left side
                    batch[k] = torch.tensor([(max_seq_len - len(f[k])) * [self.seq_pad_token_id] + f[k] if len(f[k]) < max_seq_len else f[k][-max_seq_len:] 
                                            for f in examples])
                    #batch[k] = torch.tensor([f[k] for f in examples])

                elif k == "cell_embedding":
                    # we just simply aggregate the cell embedding into fixed number of cells
                    # t_start = time.time()
                    batch[k] = torch.stack(
                        [
                            (F.adaptive_avg_pool1d(torch.tensor(f[k]).transpose(1, 0),  min_cell_cnt)).transpose(1, 0)
                            for f in examples
                        ]
                    )
                    # t_end = time.time()
                    # print(f"Collate time: {t_end - t_start}")

                else: # peak_value
                    batch[k] = torch.tensor([f[k] for f in examples])

        batch['seq_len'] = torch.tensor(_seq_len)
        return batch

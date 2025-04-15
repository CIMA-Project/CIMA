#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
@File    :		_scrna_data_processor_base.py
@Time    :   	2024/05/26 11:32:01
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
from typing import Dict, List, Optional, Union
import scanpy as sc

logger = logging.getLogger(__name__)

class SCRNAInputDataProcessorBase:
    def __init__(
        self,
        raw_adata_file_name: str,
        support_format: Optional[str]="h5ad",
        is_count_raw_data: Optional[bool]=False,
        batch_key: Optional[str]=None,
    ):
        if not os.path.isfile(raw_adata_file_name):
            raise FileNotFoundError(f"raw input scRNA .h5ad file {raw_adata_file_name} not found")
        if os.path.basename(raw_adata_file_name).split(".")[1] != support_format:
            raise ValueError(f"only {support_format} format support, but not {raw_adata_file_name} ")
        
        logger.info(f"loading {raw_adata_file_name}")
        self.adata = sc.read_h5ad(raw_adata_file_name)

        self.data_config = dict(
            input_data_file_name=raw_adata_file_name
        )

        # adata.X is raw count, otherwise raw count is in adata.raw.X
        self.data_config["is_count_raw_data"] = is_count_raw_data
        self.data_config["batch_key"] = batch_key
    
    def preprocess_data(
        self,
        gene_col: str = "gene_name",
        **kwargs
    ):
        # adapted from
        # https://github.com/microsoft/zero-shot-scfoundation/blob/main/sc_foundation_evals/data.py

        if gene_col not in self.adata.var.columns:
            logger.warning(f"{gene_col} not found in var columns, Using var_names instead")
            self.adata.var[gene_col] = self.adata.var_names.to_list()
        
        self.gene_col = gene_col
        self.data_config["gene_col"] = gene_col

        # note raw data shape
        self.data_config["raw_input__n_cells"] = self.adata.shape[0]
        self.data_config["raw_input__n_genes"] = self.adata.shape[1]

        return self._preprocess_data(**kwargs)
    
    def _preprocess_data(
        self,
        **kwargs
    ):
        raise NotImplementedError
    
    def get_scdata_config(
        self
    ):
        return self.data_config
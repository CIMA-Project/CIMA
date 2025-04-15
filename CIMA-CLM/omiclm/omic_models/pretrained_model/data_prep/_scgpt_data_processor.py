#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
@File    :		_scgpt_data_processor.py
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
import numpy as np
import scanpy as sc
import json

from scgpt.preprocess import Preprocessor

from ._scrna_data_processor_base import SCRNAInputDataProcessorBase

logger = logging.getLogger(__name__)

class ScGPTDataProcessor(SCRNAInputDataProcessorBase):
    """Adata Processor by scGPT

    The processor takes adata file with raw RNA count metrix, saved in adata.raw. 
    
    NOTE:
    You can provided HVG adata with raw RNA count metrix. In this case, you need to set `n_hvg=0` when calling `preprocess_data`.

    The  raw RNA count metrix can be read from `adata.X`, 
        or read from `adata.raw.X` if `use_raw=True` when calling `preprocess_data`.
    
    The process including:
        - filtering adata.var_names (ie, gene names, not ensembl id) based on the `gene_vocab`. 
            The `gene_vocab` is provided by scGPT. (filtering based on gene names)
        - filtering genes and cells based on counts,
        - normalizing with `target_sum=10000`,
        - log1p (optional, True or False),
        - HVG (optional, True or False),
        - binning,

    
        Parameters
        ----------
        raw_adata_file_path : str
            _description_
        is_count_raw_data : Optional[bool], optional
            _description_, by default True
        support_format : Optional[str], optional
            _description_, by default "h5ad"
    """

    def __init__(
        self,
        raw_adata_file_name: str,
        support_format: Optional[str]="h5ad",
        is_count_raw_data: Optional[bool]=False,
        **kwargs,
    ):
        
        self.is_count_raw_data = is_count_raw_data
        self.raw_adata_file_name = raw_adata_file_name
    
        super().__init__(
            raw_adata_file_name=raw_adata_file_name,
            support_format=support_format,
            is_count_raw_data=is_count_raw_data,
            **kwargs,
        )
    
    def _preprocess_data(
        self,
        gene_vocab: Union[str, List[str]],
        output_dir: Optional[str]=None,
        use_raw: bool=True,
        min_gene_vocab_matched_frac: Optional[float]=0.6,
        filter_gene_by_counts: int = 3,
        filter_cell_by_counts: Union[int, bool] = False,
        normalize_total: int = 1e4,
        n_hvg: Union[int, bool] = 1200,
        n_bins: int = 51,
        normed_key: str = "X_normed",
        log1p_key: str = "X_log1p",
        binned_key: str = "X_binned",
        **kwargs
    ):
        """_summary_

        Parameters
        ----------
        gene_vocab : Union[str, List[str]]
            _description_
        output_dir : Optional[str], optional
            _description_, by default None
        use_raw : bool, optional
            _description_, by default True
        min_gene_vocab_matched_frac : Optional[float], optional
            _description_, by default 0.6
        filter_gene_by_counts : int, optional
            _description_, by default 3
        filter_cell_by_counts : Union[int, bool], optional
            _description_, by default False
        normalize_total : int, optional
            _description_, by default 1e4
        n_hvg : Union[int, bool], optional
            NOTE: if you already have HVG, you should set `n_hvg=0`, by default 1200
        n_bins : int, optional
            _description_, by default 51
        normed_key : str, optional
            _description_, by default "X_normed"
        log1p_key : str, optional
            _description_, by default "X_log1p"
        binned_key : str, optional
            _description_, by default "X_binned"

        Raises
        ------
        ValueError
            _description_
        FileNotFoundError
            _description_
        ValueError
            _description_
        """
        if output_dir is None:
            raise ValueError("output_path for process scRNA data with Geneformer not specified")
        # Ensure save_directory exists
        if not os.path.exists(output_dir):
            # multiple GPUs will raise FileExistsError if exist_ok=False(default value)
            os.makedirs(output_dir, exist_ok=True)
        
        if isinstance(gene_vocab, str) and not os.path.isfile(gene_vocab):
            raise FileNotFoundError(f"Gene vocab for scGPT not found, {gene_vocab}")
        
        if isinstance(gene_vocab, str) and os.path.isfile(gene_vocab):
            with open(gene_vocab, "r") as f:
                _vocab_dict = json.load(f)
                _gene_vocab = list(_vocab_dict.keys())
        if isinstance(gene_vocab, (list, tuple)):
             _gene_vocab = gene_vocab

        # preprocess the data
        self.adata.var["id_in_vocab"] = [
            1 if gene in _gene_vocab else -1 
            for gene in self.adata.var[self.gene_col]
        ]

        gene_ids_in_vocab = np.array(self.adata.var["id_in_vocab"])
        fract = np.sum(gene_ids_in_vocab >= 0)/len(gene_ids_in_vocab)

        if fract < min_gene_vocab_matched_frac:
            raise ValueError(f"Only {fract*100:.2f}% genes in the dataset are in the vocabulary!")
        
        self.adata = self.adata[:, self.adata.var["id_in_vocab"] >= 0]
        self.data_config["fract_genes_in_vocab"] = fract

        # reset use_raw parameter
        if self.is_count_raw_data and use_raw:
            logger.warning(f"You set is_count_raw_data=true, but use_raw={use_raw}, "
                           f"Geneform will use adata.X as raw count matrix."
            )
            use_raw = False
        
        if use_raw:
            self.adata.X = self.adata.raw[:, self.adata.var_names].X.copy()
            # to reduce the adata size
            self.adata.raw = None

        logger.info(
            f"Matched {np.sum(gene_ids_in_vocab >= 0)}/{len(gene_ids_in_vocab)}"
            f" genes in vocabulary of size {len(_gene_vocab)}."
        )

        if n_hvg < 1:
            n_hvg = False
        # append preprocessing parameters to run config
        _prep_param = {
            "preprocesing__filter_gene_by_counts": filter_gene_by_counts,
            "preprocesing__filter_cell_by_counts": filter_cell_by_counts,
            "preprocesing__normalize_total": normalize_total,
            "preprocesing__normed_key": normed_key,
            "preprocesing__log1p_key": log1p_key,
            "preprocesing__binned_key": binned_key,
            "preprocesing__n_bins": n_bins,
            "preprocesing__n_hvg": n_hvg,
        }
        self.data_config.update(_prep_param)

        logger.info("Preprocessing data for scGPT...")
        # adapted from 
        # https://github.com/microsoft/zero-shot-scfoundation/blob/main/sc_foundation_evals/data.py#L161
        preprocessor = Preprocessor(
            # the key in adata.layers to use as raw data
            use_key = 'X',  
            # step 1
            filter_gene_by_counts = filter_gene_by_counts, 
            # step 2
            filter_cell_by_counts = filter_cell_by_counts, 
            # 3. whether to normalize the raw data and to what sum
            normalize_total = normalize_total,  
            # the key in adata.layers to store the normalized data
            result_normed_key = normed_key, 
            # 4. whether to log1p the normalized data
            log1p = True,  
            result_log1p_key = log1p_key,
            # 5. whether to subset the raw data to highly variable genes
            subset_hvg = n_hvg,  
            hvg_flavor = ("seurat_v3" 
                          if self.data_config["is_count_raw_data"] 
                          else "cell_ranger"),
            # 6. whether to bin the raw data and to what number of bins
            binning = n_bins, 
            # the key in adata.layers to store the binned data
            result_binned_key = binned_key,  
        )

        preprocessor(self.adata, batch_key = self.data_config["batch_key"])

        base_file_name = os.path.basename(self.raw_adata_file_name).split(".")[0]

        self.adata.write_h5ad(os.path.join(output_dir, 
                                               f"{base_file_name}_preprocessed.h5ad"))

        self.data_config["preprocessed__n_cells"] = self.adata.shape[0]
        self.data_config["preprocessed__n_genes"] = self.adata.shape[1]

        with open(os.path.join(output_dir, f"{base_file_name}_preprocessed_config.json"), 'w') as f:
            f.write(json.dumps(self.data_config, indent=2, sort_keys=True) + "\n")
        
        return self.adata


        
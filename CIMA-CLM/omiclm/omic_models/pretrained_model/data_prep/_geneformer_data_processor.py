#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
@File    :		_geneformer_data_processor.py
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
import pickle
import scanpy as sc
import json

from ._scrna_data_processor_base import SCRNAInputDataProcessorBase

logger = logging.getLogger(__name__)

class GeneformerDataProcessor(SCRNAInputDataProcessorBase):
    """Adata Processor by Geneformer

    The processor takes adata file with raw RNA count metrix, saved in adata.raw. NOTE: The input adata should have been processed with HVG.

    The  raw RNA count metrix can be read from adata.X, 
        or adata.raw.X if `use_raw=True` when calling `preprocess_data`.
    
    The process including:
        - sc.pp.calculate_qc_metrics() with `log1p = False`, (i.e., without log1p process),
        - filtering out genes and cells based on counts,
        - mapping gene names into ensembl ids, default using "gene_name_id_dict.pkl" provided by Geneformer. 
            Discarding unmapped genes.
    
    NOTE:
    If you provided a ready processed adata, make sure the adata has inlcuded the following preprocess steps and arttributes:
        - sc.pp.calculate_qc_metrics() with `log1p = False`.
        - filterring out genes and cells based on counts,
        - mapping gene names into ensembl id,
            NOTE: You do not have to do normalization, 
            as normalization (with target_sum=10000) will be done by `geneformer.tokenizer.TranscriptomeTokenizer.tokenize_anndata`.
        - must have adata.obs['n_counts'], which is equal to adata.obs['total_counts'] after you processed adata.
        - if you want to futher filter out cells, add adata.obs["filter_pass"] and set 1 if you want to preserve the cells,
            otherwise, set 0. In this context, only cells with adata.obs["filter_pass"]=1 will be tokenized.
            see: `geneformer.tokenizer.TranscriptomeTokenizer.tokenize_anndata`.
    
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
        is_count_raw_data: Optional[bool]=False,
        support_format: Optional[str]="h5ad",
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
        gene_name_id: Optional[Union[str, Dict[str, str]]]="gene_name_id_dict.pkl",
        has_gene_name_id_maped: bool = False,
        use_raw: Optional[bool] = True,
        output_dir: Optional[str]=None,
        min_gene_vocab_matched_frac: Optional[float]=0.6,
        filter_cell_by_genes:Optional[int]=10,
        filter_gene_by_cells:Optional[int]=10,
        **kwargs
    ):
        """_summary_

        Parameters
        ----------
        gene_name_id : Optional[Union[str, Dict[str, str]]], optional
            dictionary with gene name and esembl id mapping, by default "gene_name_id_dict.pkl".

            gene_name_id can be .json or .pkl file with gene name to ensembl id mapping, or a dict object.
            
            default "gene_name_id_dict.pkl" comes from geneformer huggingface repo:
            https://huggingface.co/ctheodoris/Geneformer/tree/main/geneformer

            NOTE:
            There are total 40248 entries in the  "gene_name_id_dict.pkl". Entries like:
            {
                'MT-TF': 'ENSG00000210049',
                'MT-RNR1': 'ENSG00000211459',
                'MT-TV': 'ENSG00000210077',
                'MT-RNR2': 'ENSG00000210082',
                'MT-TL1': 'ENSG00000209082',
                ...
                'KRT33B': 'ENSG00000131738',
                'RPS18': 'ENSG00000231500',
                'LILRA1': 'ENSG00000104974',
                'HLA-DOB': 'ENSG00000241106',
                ...
            }

        use_raw : Optional[bool], optional
            if True, scRNA raw counts metrix X can be read from adata.raw.X, by default True

        output_dir : Optional[str], optional
            directory saving the processed adata, by default None

        min_gene_vocab_matched_frac : Optional[float], optional
            min fraction of genes in the adata that can be found in vocab/gene_name_id_dict , by default 0.6

        filter_cell_by_genes : Optional[int], optional
            min number of genes a cell has, by default 10
        filter_gene_by_cells : Optional[int], optional
            min number of cells a gene has, by default 10

        Raises
        ------
        ValueError
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
        
        if not has_gene_name_id_maped:
            if isinstance(gene_name_id, dict):
                gene_name_id_map = gene_name_id
            else:
                assert os.path.isfile(gene_name_id), f"{gene_name_id} not existing."

                gene_name_id_file_ext =  os.path.basename(gene_name_id).split(".")[-1]
                if gene_name_id_file_ext == 'json':
                    with open(gene_name_id, 'r') as f:
                        gene_name_id_map = json.load(f)
                elif gene_name_id_file_ext == 'pkl':
                    with open(gene_name_id, "rb") as f:
                        gene_name_id_map = pickle.load(f)
                else:
                    raise ValueError(f"Gene name id mapp data shoulb be dict, or a file with .json or .pkl, {gene_name_id} not supported.")
            self.adata.var['ensembl_id'] = self.adata.var[self.gene_col].map(gene_name_id_map)
        else: 
            if 'ensembl_id' not in self.adata.var.columns:
                raise ValueError("adata.var does not have `ensembl_id` attribute.")
            logger.warning(f"You provided adata with ensembl id, ignoring gene name mapping.")
        
        
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
        
        # QC
        sc.pp.calculate_qc_metrics(self.adata, 
                                   percent_top=None, 
                                   log1p=False, 
                                   inplace=True,
                                   use_raw=False,
                                   layer=None)
        
        self.adata.obs['n_counts'] = self.adata.obs['total_counts']
        sc.pp.filter_cells(self.adata, min_genes=int(filter_cell_by_genes))
        sc.pp.filter_genes(self.adata, min_cells=int(filter_gene_by_cells))
        
        self.adata.var['has_ensembl_match'] = self.adata.var['ensembl_id'].notnull()

        n_all_genes = len(self.adata.var_names)
        n_matched = self.adata.var.has_ensembl_match.sum()
        fract = n_matched / n_all_genes

        if fract < min_gene_vocab_matched_frac:
            raise ValueError(f"Only {fract*100:.2f}% genes in the dataset are in the vocabulary!")
        
        base_file_name = os.path.basename(self.raw_adata_file_name).split(".")[0]
        # save the adata.var dataframe
        self.adata.var.to_csv(os.path.join(output_dir, f"{base_file_name}_var.csv"), index = False)
        
        # filter out genes that don't have a match
        self.adata = self.adata[:, self.adata.var.has_ensembl_match]
        
        # additionally, add the order of the samples, since they will be sorted
        # to speed up forward pass
        self.adata.obs['adata_order'] = self.adata.obs_names.to_list()
        self.data_config["fract_genes_in_vocab"] = fract

        if not has_gene_name_id_maped:
            logger.info(
                f"Matched {fract*100:.2f}% genes ({n_matched}/{n_all_genes})"
                f" genes in vocabulary of size {len(gene_name_id_map)}, after cell and gene filtering."
            )
        else: 
             logger.warning(
                f"Matched {fract*100:.2f}% genes ({n_matched}/{n_all_genes})"
                f" genes after cell and gene filtering."
                f" Note, you provided adata with ensembl id, which might be potential issue"
                f" as the provided ensembl ids would be different from ensembl ids in the vocab provided by Genformer."
            )

        self.adata.write_h5ad(os.path.join(output_dir, 
                                               f"{base_file_name}_preprocessed.h5ad"))

        self.data_config["preprocessed__n_cells"] = self.adata.shape[0]
        self.data_config["preprocessed__n_genes"] = self.adata.shape[1]

        with open(os.path.join(output_dir, f"{base_file_name}_preprocessed_config.json"), 'w') as f:
            f.write(json.dumps(self.data_config, indent=2, sort_keys=True) + "\n")
        
        return self.adata
        
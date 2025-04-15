#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
@File    :		_geneformer_tokenizer.py
@Time    :   	2024/05/26 15:50:22
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
import pickle
from shutil import copyfile
from typing import Dict, List, Optional, Tuple
import json

from transformers.tokenization_utils import PreTrainedTokenizer, AddedToken
from datasets import load_from_disk

from geneformer.tokenizer import TranscriptomeTokenizer

VOCAB_FILES_NAMES = {
    "vocab_file": "token_dictionary.pkl",
    "gene_name_id_file": "gene_name_id_dict.pkl",
    "gene_median_file": "gene_median_dictionary.pkl",
}

logger = logging.getLogger(__name__)

class GeneformerTokenizerWrapper(PreTrainedTokenizer):
    vocab_files_names = VOCAB_FILES_NAMES

    def __init__(
        self,
        vocab_file: str,
        gene_name_id_file: str = None,
        gene_median_file: str = None,
        pad_token: str = "<pad>",
        mask_token: str = "<mask>",
        **kwargs,
    ):
        """_summary_

        Parameters
        ----------
        vocab_file : str
            gene token vocab file
        gene_name_id_file : str, optional
            gene names to ensembl ids mapping file, by default None
        gene_median_file : str, optional
            Path to pickle file containing dictionary of non-zero median
            gene expression values across Genecorpus-30M. by default None
        pad_token : str, optional
            by default "<pad>"
        mask_token : str, optional
            by default "<mask>"
        """
        self.vocab_file = vocab_file
        self.gene_name_id_file = gene_name_id_file
        self.gene_median_file = gene_median_file

        pad_token = AddedToken(pad_token, lstrip=False, rstrip=False) if isinstance(pad_token, str) else pad_token
        mask_token = AddedToken(mask_token, lstrip=False, rstrip=False) if isinstance(mask_token, str) else mask_token
        
        with open(vocab_file, "rb") as f:
            self.vocab = pickle.load(f)

        logger.info(f"vocab size = {len(self.vocab)}")
        self._id2vocab={v:k for k, v in self.vocab.items()}

        super().__init__(
            add_prefix_space=False,
            pad_token=pad_token,
            mask_token=mask_token,
            **kwargs,
        )
    
    #
    def __call__(
        self,
        adata_full_file_name: str,
        save_dataset_dir: str,
        cell_type_col: str = "cell_type",
        columns_to_keep: List[str] = ["adata_order"],
        num_workers: int = 1,
        **kwargs
    ):
        # get the extension from adata_path
        # _, ext = os.path.splitext(adata_full_file_name)
        # ext = ext.strip(".")
        adata_fanme = os.path.basename(adata_full_file_name).split(".")
        adata_name = adata_fanme[0]
        ext = adata_fanme[-1]

        if ext not in ["loom", "h5ad"]:
            raise ValueError(f"adata_path must be a loom or h5ad file. Got {ext}")
        if ext == "h5ad":
            msg = ("using h5ad file. This sometimes causes issues. "
                   "If not working try with loom.")
            logger.warning(msg)
        
        cols_to_keep = dict(zip([cell_type_col] + columns_to_keep, 
                                [cell_type_col] + columns_to_keep))
        

        # initialize tokenizer
        geneformor_tokenizer = TranscriptomeTokenizer(
            cols_to_keep, 
            nproc = num_workers,
            # gene_median_file: Path to pickle file containing dictionary of non-zero median
            # gene expression values across Genecorpus-30M.
            gene_median_file=self.gene_median_file,
            token_dictionary_file=self.vocab_file,
            model_input_size=2048,  # default value
            special_token=False,    # default value
        )


        # get the top directory of the adata_path
        adata_dir = os.path.dirname(adata_full_file_name)
        geneformor_tokenizer.tokenize_data(
            adata_dir,
            save_dataset_dir, 
            adata_name,
            file_format=ext
        )

        return load_from_disk(os.path.join(save_dataset_dir, f"{adata_name}.dataset"))

    @property
    def vocab_size(self):
        return len(self.vocab)
    
    def get_vocab(self) -> Dict[str, int]:
        return self.vocab
    
    def _convert_token_to_id(self, token):
        return self.vocab[token]
       
    def _convert_id_to_token(self, index:int):
        return self._id2vocab[index]
    
    def save_vocabulary(self, save_directory: str, filename_prefix: Optional[str] = None) -> Tuple[str]:
        if not self.can_save_slow_tokenizer:
            raise ValueError(
                "Your fast tokenizer does not have the necessary information to save the vocabulary for a slow "
                "tokenizer."
            )

        if not os.path.isdir(save_directory):
            logger.error(f"Vocabulary path ({save_directory}) should be a directory")
            return
        
        out_vocab_file = os.path.join(
            save_directory, (filename_prefix + "-" if filename_prefix else "") + VOCAB_FILES_NAMES["vocab_file"]
        )

        out_gene_name_id_file = os.path.join(
            save_directory, (filename_prefix + "-" if filename_prefix else "") + VOCAB_FILES_NAMES["gene_name_id_file"]
        )

        out_gene_median_file = os.path.join(
            save_directory, (filename_prefix + "-" if filename_prefix else "") + VOCAB_FILES_NAMES["out_gene_median_file"]
        )

        if os.path.abspath(self.vocab_file) != os.path.abspath(out_vocab_file):
            copyfile(self.vocab_file, out_vocab_file)
            with open(os.path.splitext(out_vocab_file) + '.json', 'w') as f:
                f.write(json.dumps(self.vocab, indent=2) + "\n")

        if os.path.abspath(self.gene_name_id_file) != os.path.abspath(out_gene_name_id_file):
            copyfile(self.gene_name_id_file, out_gene_name_id_file)
            with open(self.gene_name_id_file, 'rb') as f:
                gene_name_id = pickle.load(f)
            with open(os.path.splitext(out_gene_name_id_file) + '.json', 'w') as f:
                f.write(json.dumps(gene_name_id, indent=2) + "\n")
        
        if os.path.abspath(self.gene_name_id_file) != os.path.abspath(out_gene_name_id_file):
            copyfile(self.gene_median_file, out_gene_median_file)
            with open(self.gene_median_file, 'rb') as f:
                gene_median = pickle.load(f)
            with open(os.path.splitext(out_gene_median_file) + '.json', 'w') as f:
                f.write(json.dumps(gene_median, indent=2) + "\n")

        return (out_vocab_file, out_gene_name_id_file, out_gene_median_file)
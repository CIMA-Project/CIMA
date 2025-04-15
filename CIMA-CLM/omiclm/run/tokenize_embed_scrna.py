#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
@File    :		tokenize_embed_scrna.py
@Time    :   	2024/05/24 13:08:37
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
import sys
from typing import Any, Dict, List, Mapping, NewType
import logging
from pathlib import Path
import transformers
import datasets
from dataclasses import dataclass
import scanpy as sc

from transformers import (
    HfArgumentParser,
    TrainingArguments,
    Trainer,
    set_seed
)

from datasets import Dataset, DatasetDict, load_dataset, load_from_disk
from transformers.utils import check_min_version

sys.path.append(
    str(Path(__file__).resolve().parents[1]) #
)
from omic_config import (
    OmicRawDataConfig,
    OmicPretrainedModelAndTokenizationConfig,
    OmicFormerTrainingArguments,
)

from omic_models import (
    PretrainedModelName,
    PRETRAINED_MODEL_NAME_CLS_MAP,
    PRETRAINED_TOKENIZER_NAME_CLS_MAP,
    SCRNA_DATA_PREPROCESSOR_FOR_PRETRAINED_MODEL_MAP,
)

check_min_version("4.39.3")
logger = logging.getLogger(__name__)


def main():
    parser = HfArgumentParser((
        OmicRawDataConfig, 
        OmicPretrainedModelAndTokenizationConfig, 
        OmicFormerTrainingArguments))
    raw_data_config, pretrained_model_tokenizer_config, training_config = parser.parse_args_into_dataclasses()

    # Setup logging
    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(name)s - %(message)s",
        datefmt="%m/%d/%Y %H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout)],
    )

    if training_config.should_log:
        # The default of training_args.log_level is passive, 
        # so we set log level at info here to have that default.
        transformers.utils.logging.set_verbosity_info()

    log_level = training_config.get_process_log_level() # log_level = 20
    logger.setLevel(log_level)
    datasets.utils.logging.set_verbosity(log_level)
    transformers.utils.logging.set_verbosity(log_level)
    transformers.utils.logging.enable_default_handler()
    transformers.utils.logging.enable_explicit_format()

    # Log on each process the small summary:
    logger.warning(
        f"Process rank: {training_config.local_rank}, device: {training_config.device}, n_gpu: {training_config.n_gpu}, "
        + f"distributed training: {training_config.parallel_mode.value == 'distributed'}, 16-bits training: {training_config.fp16}"
    )

    #################################################################
    # loading row dataset and preprocess
    #################################################################

    raw_scrna_count_file_name = raw_data_config.raw_scrna_file_name
    sc_model_dir = pretrained_model_tokenizer_config.scrna_model_path
    logger.info(f"sc model pretrained path:{sc_model_dir}")

    logger.info(f"loading scRNA data from {raw_scrna_count_file_name}") 
    processed_cell_type = raw_data_config.processed_cell_type_name

    adata = sc.read_h5ad(raw_scrna_count_file_name)
    if 'log1p' in adata.uns_keys():
        logger.info("remove 'log1p' key from the adata.uns")
        adata.uns.pop('log1p')

    logger.info(f"raw scRNA data:\n{adata}")

    # preprocess scRNA dataset
    adata_preprocessor = SCRNA_DATA_PREPROCESSOR_FOR_PRETRAINED_MODEL_MAP[pretrained_model_tokenizer_config.scrna_model_name](
        raw_adata_file_name=raw_scrna_count_file_name
    )
    prep_out_dir = os.path.join(raw_data_config.raw_scrna_preprocessed_dir, processed_cell_type)
    
    scrna_preprocess_args = dict(
        gene_name_id=os.path.join(pretrained_model_tokenizer_config.scrna_model_path, "gene_name_id_dict.pkl"),
        has_gene_name_id_maped=False,
        use_raw=True,
        output_dir=prep_out_dir,
        min_gene_vocab_matched_frac=0.6,
        filter_cell_by_genes=10,
        filter_gene_by_cells=10,
        # <-- geneforerm preprocess args
        # >-- scgpt preprocess args (removed duplicated args)
        gene_vocab=os.path.join(sc_model_dir, "vocab.json"),
        filter_gene_by_counts=3,
        filter_cell_by_counts=False,
        normalize_total=1e4,
        n_hvg=0 if raw_data_config.raw_scrna_has_hvg else 5000,
        n_bins=51,
        normed_key="X_normed",
        log1p_key="X_log1p",
        binned_key="X_binned",
    )

    # scRNA data preprocess
    logger.info(f"preprocessing scRNA data with {pretrained_model_tokenizer_config.scrna_model_name.value}")
    adata_processed = adata_preprocessor.preprocess_data(**scrna_preprocess_args)

    #################################################################
    # Tokenize (if neccessary) and embedding 
    #################################################################
    # Geneformer need tokenization before embedding
    tokenized_dataset = None
    if pretrained_model_tokenizer_config.scrna_model_name == PretrainedModelName.GENEFORMER:
        with training_config.main_process_first(desc="scRNA data tokenization by Geneformer"):
            logger.info(f"scRNA tokenization...")
            tokenizer = PRETRAINED_TOKENIZER_NAME_CLS_MAP[pretrained_model_tokenizer_config.scrna_model_name].from_pretrained(
                sc_model_dir, 
                local_files_only=True)
            
            tokenized_dataset = tokenizer(
                adata_full_file_name = os.path.join(prep_out_dir, f"{processed_cell_type}_preprocessed.h5ad"),
                save_dataset_dir=pretrained_model_tokenizer_config.tokenized_scrna_dataset_dir,
                cell_type_col=raw_data_config.raw_scrna_cell_type_var_name,
                columns_to_keep=["sample", "adata_order"],
                num_workers = 4,
            )
    
    logger.info(f"scRNA embedding...")
    sc_model_cls =  PRETRAINED_MODEL_NAME_CLS_MAP[pretrained_model_tokenizer_config.scrna_model_name]

    sc_model_init_args = dict(
        pretrained_model_name_or_path=sc_model_dir,
        output_attentions=False,
        output_hidden_states=True,
        local_files_only=True,
        # <-- geneforerm model args
        # >-- scgpt model args (removed duplicated args)
        config_file_name="args.json",
        vocab_file_name="vocab.json", 
        model_file_name="best_model.pt",

        device = training_config.device
    )

    sc_model_embedding_args = dict(
        tokenized_dataset=tokenized_dataset,
        adata=adata_processed,
        batch_size = 48,
        embedding_key = "cell_embedding",
        layer=-2,
        pad_token_id=0,
        # <-- geneforerm embedding args
        # >-- scgpt embedding args (removed duplicated args)
        adata_or_file=adata_processed,
        gene_col="gene_name",
        max_length=5000, # for scGPT
        cell_embedding_mode="cls", 
        obs_to_save = "sample",
        return_new_adata=True,  
    )

    sc_model = sc_model_cls.from_pretrained(**sc_model_init_args)
    embeded_adata = sc_model.extract_sample_embedding(**sc_model_embedding_args)

    save_dir = os.path.join(pretrained_model_tokenizer_config.emb_scrna_dataset_dir, processed_cell_type)
    logger.info(f"save scRNA embedding at {save_dir}")

    embeded_adata.write(save_dir + f"/{processed_cell_type}_embedding.h5ad")

    # embedding = embeded_adata.obsm['cell_embedding'].tolist()
    # sample = list(embeded_adata.obs['sample'])
    # cell_id = (list(embeded_adata.obs['adata_order']) 
    #            if pretrained_model_tokenizer_config.scrna_model_name == PretrainedModelName.GENEFORMER 
    #            else list(embeded_adata.obs_names))
    
    # output_ds = Dataset.from_dict(
    #     {
    #         'sample': sample,
    #         'cell_id': cell_id,
    #         'cell_embedding': embedding,
    #     }
    # )
    # output_ds.save_to_disk(save_dir)
    # logger.info("DONE")


if __name__=='__main__':
    main()


    
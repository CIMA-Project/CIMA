#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
@File    :		run_omicformer_prediction.py
@Time    :   	2024/07/03 16:59:03
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

# os.environ["CUDA_VISIBLE_DEVICES"] = "2,3"

import sys
from typing import Any, Dict, List, Mapping, NewType
import logging
from pathlib import Path
import json
import math
import copy
import gc

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
import torch.utils.data
import transformers
import datasets
import evaluate
from tqdm import tqdm
from dataclasses import dataclass
import scanpy as sc

from sklearn.metrics import mean_squared_error
from scipy import stats

import time

from transformers import (
    HfArgumentParser,
    EvalPrediction,
    set_seed,
    TrainingArguments
)

from datasets import Dataset, DatasetDict, load_dataset, load_from_disk
from transformers.utils import check_min_version
from transformers.trainer_utils import get_last_checkpoint

sys.path.append(
    str(Path(__file__).resolve().parents[1]) #
)

from omic_config import (
    OmicRawDataConfig,
    OmicPretrainedModelAndTokenizationConfig,
    OmicFormerTrainingArguments,
    OmicFormerTrainer,
)

from omic_models import (
    OmicFormerConfig, 
    OmicFormerPreTrainedModel,

    PretrainedModelName,
    PRETRAINED_MODEL_NAME_CLS_MAP,
    PRETRAINED_TOKENIZER_NAME_CLS_MAP,
)

from omic_dataset_collator import OmicDataset, OmicDataCollator

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

    set_seed(training_config.seed)

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
    # load TOKENIZED SEQ dataset for prediction
    #################################################################    
    with training_config.main_process_first(desc="loading seq data and model"):

        seq_tokenized_dataset_dir = pretrained_model_tokenizer_config.tokenized_seq_dataset_dir + f"/{raw_data_config.processed_cell_type_name}"
        logger.info(f">>> loading TOKENIZED SEQ data and split from {seq_tokenized_dataset_dir}")
        
        seq_hf_ds = load_from_disk(seq_tokenized_dataset_dir)

        # this is only for precidtion
        tst_dataset = seq_hf_ds['test']
        # tst_dataset =  tst_dataset.select(range(1000))
        logger.info(f"test dataset: \n{tst_dataset}")
        # Dataset({
        #     features: ['sample', 'pos', 'peak_value', 'sum_peak_value', 'record_id', 'input_ids', 'norm_peak_value', 'log1p_norm_peak_value'],
        #     num_rows: 6756697
        # })

    scrna_embedding_dataset_dir = os.path.join(pretrained_model_tokenizer_config.emb_scrna_dataset_dir, raw_data_config.processed_cell_type_name)
    
    # logger.info(f">>> SCRNA EMBEDDING dataset:\n{scrna_embedding_dataset_trn}")
    ac_data_file_name = os.path.join(scrna_embedding_dataset_dir, f"{raw_data_config.processed_cell_type_name}_embedding.h5ad")
    logger.info(f">>> loading SCRNA EMBEDDING data (h5ad) from {ac_data_file_name}")
    sc_adata = sc.read_h5ad(ac_data_file_name)
    
    tst_pt_ds: torch.utils.data.Dataset = OmicDataset(
        tst_dataset, 
        sc_adata,
        peak_value_feature=pretrained_model_tokenizer_config.peak_value_col_name,
    )

    ####################################################################
    # SEQ embedding model
    #
    ####################################################################
    logger.info(f">>> LOADED pretrained model and tokenizer from {pretrained_model_tokenizer_config.seq_model_path}")
    seq_emb_model = (PRETRAINED_MODEL_NAME_CLS_MAP[
        pretrained_model_tokenizer_config.seq_model_name
        ].from_pretrained(pretrained_model_tokenizer_config.seq_model_path, use_cache=True)
    ) if not pretrained_model_tokenizer_config.seq_input_has_embedding else None
    seq_emb_model.requires_grad_(False)

    # logger.info(f">>> LOADED pretrained model and tokenizer from {pretrained_model_tokenizer_config.scrna_model_path}")
    # scrna_emb_model = (PRETRAINED_MODEL_NAME_CLS_MAP[
    #     pretrained_model_tokenizer_config.scrna_model_name
    #     ].from_pretrained(pretrained_model_tokenizer_config.scrna_model_path, use_cache=True)
    # ) if not pretrained_model_tokenizer_config.scrna_input_has_embedding else None
    
    seq_input_pooling_size = 501    # max seq length for seq tokenized data of `input_ids`
    scrna_input_pooling_size = 200  # pooling size of cell conts after scrna embeddings
    
    omic_data_collator = OmicDataCollator(
        seq_pad_token_id=4, 
        sc_max_cnt=scrna_input_pooling_size, 
        seq_max_len=seq_input_pooling_size,
        seq_emb_model=seq_emb_model,
    )
    
    omic_model = OmicFormerPreTrainedModel.from_pretrained(
        pretrained_model_name_or_path=pretrained_model_tokenizer_config.omicformer_pretrained_model_dir,
        local_files_only=True,
        
        seq_emb_model=seq_emb_model,
        seq_emb_dim=pretrained_model_tokenizer_config.seq_model_d_dim,
        seq_emb_extraction_kwargs=None,

        scrna_emb_model=None,
        scrna_emb_dim=pretrained_model_tokenizer_config.scrna_model_d_dim, 
        scrna_emb_extraction_kwargs=None,
    ).to(training_config.device)
    
    omic_model.config.use_cache = False
    logger.info(f">>> OmicFormerPreTrainedModel: \n{omic_model}")

    logger.info(f"num params: {omic_model.num_parameters()}")
    logger.info(f"num trainable params: {omic_model.num_parameters(only_trainable=True)}")

    # pred_args = TrainingArguments(
    #     output_dir = training_config.output_dir + '/prediction',
    #     do_train = False,
    #     do_predict = True,
    #     per_device_eval_batch_size = training_config.per_device_eval_batch_size,    
    #     dataloader_drop_last = False,
    #     report_to=None,
    #     dataloader_num_workers = training_config.dataloader_num_workers,
    #     disable_tqdm=training_config.disable_tqdm,
    #     # eval_accumulation_steps=training_config.eval_accumulation_steps, # see https://discuss.huggingface.co/t/batch-size-for-trainer-predict/3374
    #     ddp_timeout=training_config.ddp_timeout,
    # )

    # trainer = OmicFormerTrainer(
    #     model=omic_model,
    #     args=pred_args,
    #     eval_dataset=tst_pt_ds,
    #     data_collator=omic_data_collator,
    # )
    # pred_out = trainer.predict(tst_pt_ds)

    data_loader = torch.utils.data.DataLoader(
        tst_pt_ds,
        batch_size=256,
        collate_fn=omic_data_collator,
        num_workers=3,
        pin_memory=True,
        drop_last=False,
    )

    omic_model.eval()
    pred_results = dict(
        sample_record_id=[],
        #peak_value=[],
        pred_peak_value=[],
    )
    for batch in data_loader:

        seq_input_ids = batch['seq_input_ids'].to(training_config.device)
        seq_len = batch['seq_len'].to(training_config.device)
        cell_embedding = batch['cell_embedding'].to(training_config.device)
        with torch.no_grad():
            pred_out = omic_model(seq_input_ids=seq_input_ids, seq_len=seq_len, cell_embedding=cell_embedding)
        
        if isinstance(pred_out, tuple):
            pred_out = pred_out[0]
        elif isinstance(pred_out, Mapping):
            pred_out = pred_out['logits']

        pred_out = pred_out.squeeze().cpu().numpy().tolist()
        # peak_value = batch['peak_value'].numpy().tolist()
        # sample_record_id = batch['sample_record_id'].numpy().tolist()
        # pred_results['sample_record_id'].extend(sample_record_id)
        # pred_results['peak_value'].extend(peak_value)
        sample_record_id = batch['sample_record_id'].numpy().tolist()
        pred_results['sample_record_id'].extend(sample_record_id)
        pred_results['pred_peak_value'].extend(pred_out)

    pred_results_dir = training_config.output_dir + '/prediction'
    if not os.path.exists(pred_results_dir):
        os.makedirs(pred_results_dir, exist_ok=True)

    with training_config.main_process_first(desc="clean scRNA cache files"):
        result_save_path = os.path.join(pred_results_dir, 'silicon_results.csv')
        logger.info(f">>> saving prediction results to {result_save_path}")
        pd.DataFrame.from_dict(pred_results, orient='columns').to_csv(result_save_path, index=False)
        
        seq_hf_ds.cleanup_cache_files()

    logger.info("<<<<<<<<<<<<<<<<Done")


if __name__=='__main__':
    main()



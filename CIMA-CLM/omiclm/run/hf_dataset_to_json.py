#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
@File    :		hf_dataset_to_json.py
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
import json
import gzip
import multiprocessing

import numpy as np
import torch
import torch.nn.functional as F
import transformers
import datasets
from tqdm import tqdm
from dataclasses import dataclass

from datasets.utils import tqdm as hf_tqdm
 

from transformers import (
    HfArgumentParser,
    TrainingArguments,
    Trainer,
    DefaultDataCollator,
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
        OmicPretrainedModelAndTokenizationConfig))
    raw_data_config, pretrained_model_tokenizer_config = parser.parse_args_into_dataclasses()

    training_config = TrainingArguments(output_dir="./tmp_dir")

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
    # load TOKENIZED dataset
    #################################################################

    tokenized_dataset_dir = pretrained_model_tokenizer_config.tokenized_seq_dataset_dir + f"/{raw_data_config.processed_cell_type_name}"
    logger.info(f"loading seq tokenized dataset from: {tokenized_dataset_dir}")

    # columns_to_remove = ['seq', 'attention_mask']

    with training_config.main_process_first(desc="loading data and model"):
        seq_hf_ds = load_from_disk(tokenized_dataset_dir)

        logger.info(f"seq tokenized datasets: \n{seq_hf_ds}")
        # >>> SEQ tokenized data:
        # DatasetDict({
        #     train: Dataset({
        #         features: ['sample', 'pos', 'seq', 'peak_value', 'input_ids', 'attention_mask', 'record_id'],
        #         num_rows: 124167712
        #     })
        #     test: Dataset({
        #         features: ['sample', 'pos', 'seq', 'peak_value', 'input_ids', 'attention_mask', 'record_id'],
        #         num_rows: 4208429
        #     })
        # })

        train_val_dataset = seq_hf_ds['train'].train_test_split(test_size=0.0005, shuffle=True, seed=42)
        
        trn_dataset = train_val_dataset['train']
        logger.info(f"train dataset: \n{trn_dataset}")
        # Dataset({
        #     features: ['sample', 'pos', 'seq', 'peak_value', 'input_ids', 'attention_mask', 'record_id'],
        #     num_rows: 111750940
        # })

        val_dataset = train_val_dataset['test']
        logger.info(f"validation dataset: \n{val_dataset}")
        # Dataset({
        #     features: ['sample', 'pos', 'seq', 'peak_value', 'input_ids', 'attention_mask', 'record_id'],
        #     num_rows: 12416772
        # })

        tst_dataset = seq_hf_ds['test']
        logger.info(f"test dataset: \n{tst_dataset}")
        # Dataset({
        #     features: ['sample', 'pos', 'seq', 'peak_value', 'input_ids', 'attention_mask', 'record_id'],
        #     num_rows: 4208429
        # })

    ############################################################################
    # for debug
    ############################################################################
    # if t_args.do_predict:
    #     trn_dataset = trn_dataset.select(range(32))
    # # update value of dataset
    # specific_sample = trn_dataset[0]
    # # Modify the length of the list in the specific sample
    # specific_sample['input_ids'] = specific_sample['input_ids'][:100]
    # # Replace the sample in the dataset
    # trn_dataset = trn_dataset.map(lambda x, idx: specific_sample if idx == 0 else x, with_indices=True)
    # for ds in trn_dataset:
    #     logger.info(f"{ds['record_id']}, {len(ds['input_ids'])}")

    dstype = ['train', 'validation', 'test']
    tokenized_ds = [trn_dataset, val_dataset, tst_dataset]
    
    n_shard = 500 # i.e, the number of generated json files for each dataset type, like chunks
    for _dstype, _token_ds in zip(dstype, tokenized_ds):
        _shard_save_dir = os.path.join(f"{tokenized_dataset_dir}_json", _dstype)
        if not os.path.exists(_shard_save_dir):
            os.makedirs(_shard_save_dir, exist_ok=True)

        for i_shardx in range(0, n_shard):
            i_dataset = _token_ds.shard(num_shards=n_shard, index=i_shardx)
            i_file_name = os.path.join(_shard_save_dir, f"{_dstype}_{i_shardx}_{n_shard}.json")
            i_dataset.to_json(i_file_name, num_proc=50, batch_size=10000, orient='records', lines=True)
            logger.info(f"sharded {i_shardx}/{n_shard} from {_dstype}")
        # with multiprocessing.Pool(5) as pool:
        #     for out_str in hf_tqdm(
        #             pool.imap(
        #                 get_shard,
        #                 [(n_shard, i_shardx, _token_ds, _dstype, _shard_save_dir) for i_shardx in range(0, n_shard)],
        #             ),
        #             total=n_shard,
        #             unit="ba",
        #             desc="Writing dataset shard",
        #         ):
        #         logger.info(f"sharded {out_str} from {_dstype}")
    
    logger.info("clean up train_test_split cache files.")
    with training_config.main_process_first(desc="cleaning up cache files"):
        train_val_dataset.cleanup_cache_files()
        tst_dataset.cleanup_cache_files()
    logger.info("DONE")

def get_shard(args):
    n_shard, i_shard, ds, _dstype, _shard_save_dir = args
    i_dataset = ds.shard(num_shards=n_shard, index=i_shard)
    i_file_name = os.path.join(_shard_save_dir, f"{_dstype}_{i_shard}_{n_shard}.json")
    i_dataset.to_json(i_file_name, num_proc=10, batch_size=10000, orient='records', lines=True)
    return f"{i_shard}/{n_shard}"

if __name__=='__main__':
    main()


    



#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
@File    :		build_datasets.py
@Time    :   	2024/05/16 12:44:59
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
import argparse
import logging
from pathlib import Path
import shutil
import transformers
import datasets
import pandas as pd
import numpy as np

from transformers import (
    HfArgumentParser,
    TrainingArguments,
    set_seed
)

from datasets import load_dataset, DatasetDict, load_from_disk
from transformers.utils import check_min_version

sys.path.append(
    str(Path(__file__).resolve().parents[1]) #
)
from omic_config import OmicRawDataConfig, OmicPretrainedModelAndTokenizationConfig

check_min_version("4.39.3")
logger = logging.getLogger(__name__)

def main():

    parser = HfArgumentParser((OmicRawDataConfig, OmicPretrainedModelAndTokenizationConfig))
    raw_ds_config, pretrained_model_tokenizer_config = parser.parse_args_into_dataclasses()

    # Setup logging
    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(name)s - %(message)s",
        datefmt="%m/%d/%Y %H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout)],
    )

    transformers.utils.logging.set_verbosity_info()

    log_level = 20 # log_level = 20
    logger.setLevel(log_level)
    datasets.utils.logging.set_verbosity(log_level)
    transformers.utils.logging.set_verbosity(log_level)
    transformers.utils.logging.enable_default_handler()
    transformers.utils.logging.enable_explicit_format()

    seq_tokenized_dataset_dir = pretrained_model_tokenizer_config.tokenized_seq_dataset_dir + f"/{raw_ds_config.processed_cell_type_name}"
    
    # first copy the original tokenized SEQ data to the tmp_dir
    #tmp_dir = seq_tokenized_dataset_dir + '_original_copied'
    # tmp_dir = "pythonproject/DNALLM/datasets/tokenized_datasets/seq_dataset/hyenadna/CD8_expression_5K_original_copied"
    tmp_dir = "project/DNALLM/datasets/tokenized_datasets/seq_dataset/hyenadna/ncMono_C1QA_original_copied"
    # if os.path.isdir(tmp_dir):
    #     try:
    #         shutil.rmtree(tmp_dir, ignore_errors=True)
    #     except Exception as e:
    #         raise ValueError(f"Failed to remove: {e}")
    # os.makedirs(tmp_dir, exist_ok=True)
    
    # logger.info(f">>> copying original tokenized SEQ data from {seq_tokenized_dataset_dir} to {tmp_dir}")
    # for item in os.listdir(seq_tokenized_dataset_dir):
    #     item_path = os.path.join(seq_tokenized_dataset_dir, item)
    #     if item_path == os.path.abspath(seq_tokenized_dataset_dir):
    #         continue
    #     if os.path.isfile(item_path):
    #         shutil.copy2(item_path, tmp_dir)
    #     else:
    #         dst_dir = os.path.join(tmp_dir, item)
    #         os.makedirs(dst_dir, exist_ok=True)
    #         shutil.copytree(item_path, dst_dir, dirs_exist_ok=True)
    
    logger.info(f">>> removing original tokenized SEQ data from {seq_tokenized_dataset_dir}")
    try:
        shutil.rmtree(seq_tokenized_dataset_dir, ignore_errors=True)
        os.makedirs(seq_tokenized_dataset_dir, exist_ok=True)
    except Exception as e:
        raise ValueError(f"Failed to remove: {e}")
    
    logger.info(f">>> loading TOKENIZED SEQ data and split from {tmp_dir}")
    seq_hf_ds = load_from_disk(tmp_dir)

    peak_value_file = pretrained_model_tokenizer_config.normalized_peak_value_file

    logger.info(f">>> loading normalized peak value from {peak_value_file}")
    peak_values = pd.read_csv(peak_value_file, index_col='sample', chunksize=20)
    peak_value_df = pd.concat(peak_values, axis=0)

    def add_peak_value(example):
        #example['norm_peak_value'] = (peak_value_df.loc[example['sample'], example["pos"]]).astype(np.float32)
        #example['log1p_norm_peak_value'] = np.log1p(example['norm_peak_value']).astype(np.float32)
        # example.pop('seq', None)
        # example.pop('attention_mask', None)
        #return example
        try:
            # 尝试获取 norm_peak_value
            norm_peak_value = peak_value_df.loc[example['sample'], example["pos"]].astype(np.float32)
            example['norm_peak_value'] = norm_peak_value
            example['log1p_norm_peak_value'] = np.log1p(example['norm_peak_value']).astype(np.float32)
            return example

        except KeyError:
            #logger.info(f">>> can not loc {example['sample']} and {example['pos']}")
            norm_peak_value = np.float32(0.0)
            example['norm_peak_value'] = norm_peak_value
            example['log1p_norm_peak_value'] = np.log1p(example['norm_peak_value']).astype(np.float32)
            example['sample'] = "None"
            return example

    def filter_none_samples(example):
        return example['sample'] != "None"


    running_main_process_config = TrainingArguments(output_dir="./tmp_dir")

    with running_main_process_config.main_process_first(desc="adding peak value for training tokenized seq data"):
        trn_ds = seq_hf_ds['train'].map(
            add_peak_value,
            num_proc=30,
            remove_columns=['seq', 'attention_mask']
        )
        trn_ds = trn_ds.filter(filter_none_samples,num_proc=30)


    with running_main_process_config.main_process_first(desc="adding peak value for test tokenized seq data"):
        tst_ds = seq_hf_ds['test'].map(
            add_peak_value,
            num_proc=30,
            remove_columns=['seq', 'attention_mask']
        )
        tst_ds = tst_ds.filter(filter_none_samples,num_proc=30)

    res_ds = DatasetDict()
    res_ds['train'] = trn_ds
    res_ds['test'] = tst_ds

    logger.info(f">>> saving tokenized SEQ data with normalized peak value to {seq_tokenized_dataset_dir}")
    res_ds.save_to_disk(seq_tokenized_dataset_dir)

    with running_main_process_config.main_process_first(desc="cleaning up cache files"):
        res_ds.cleanup_cache_files()
        seq_hf_ds.cleanup_cache_files()

    logger.info(">>>Done!")

if __name__=='__main__':
    main()


    



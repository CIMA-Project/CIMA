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
import sys
import argparse
import logging
from pathlib import Path
import transformers
import datasets
from transformers import (
    HfArgumentParser,
    set_seed
)

from datasets import load_dataset, DatasetDict
from transformers.utils import check_min_version

sys.path.append(
    str(Path(__file__).resolve().parents[1]) #
)
from omic_config import OmicRawDataConfig

check_min_version("4.39.3")
logger = logging.getLogger(__name__)

def main():
    parser = HfArgumentParser((OmicRawDataConfig))
    raw_ds_config, = parser.parse_args_into_dataclasses()

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

    if raw_ds_config.raw_seq_plain_file_name_for_train is None:
        raise ValueError(f"raw data in plain files for train {raw_ds_config.raw_seq_plain_file_name_for_train} not specified.")
    
    if raw_ds_config.raw_seq_plain_file_name_for_test is None:
        raise ValueError(f"raw data in plain files for test {raw_ds_config.raw_seq_plain_file_name_for_test} not specified.")

    logger.info(f"Start generate dataset for: \n{raw_ds_config.raw_seq_plain_file_name_for_train}\n{raw_ds_config.raw_seq_plain_file_name_for_test}")
    hf_ds = load_dataset("csv", 
        delimiter="\t+|\s+|,",
        name="raw_seq_dataset",  
        trust_remote_code=True, 
        data_files = {
            "train": raw_ds_config.raw_seq_plain_file_name_for_train,
            "test": raw_ds_config.raw_seq_plain_file_name_for_test
        },
        num_proc=50,
    )
    logger.info(f"generted HF datasets: \n{hf_ds}")
    # DatasetDict({
    #     train: Dataset({
    #         features: ['sample', 'pos', 'seq', 'peak_value'],
    #         num_rows: 124167712
    #     })
    #     test: Dataset({
    #         features: ['sample', 'pos', 'seq', 'peak_value'],
    #         num_rows: 4208429
    #     })
    # })

    trn_ds = hf_ds['train']
    trn_len = len(trn_ds)
    trn_ds = trn_ds.add_column('record_id', list(range(trn_len)))
    logger.info(f"train dataset: \n{trn_ds}")

    tst_ds = hf_ds['test']
    tst_len = len(tst_ds)
    tst_ds = tst_ds.add_column('record_id', list(range(trn_len, trn_len + tst_len)))
    logger.info(f"test dataset: \n{tst_ds}")

    res_ds = DatasetDict()
    res_ds['train'] = trn_ds
    res_ds['test'] = tst_ds

    res_ds.save_to_disk(raw_ds_config.raw_seq_dataset_dir)
    res_ds.cleanup_cache_files()
    logger.info(f"HF datasets saved at: \n{raw_ds_config.raw_seq_dataset_dir}")

if __name__=='__main__':
    main()


    



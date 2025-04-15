#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
@File    :		tokenize_datasets.py
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

import sys
import logging
from pathlib import Path
import transformers
import datasets
from transformers import (
    HfArgumentParser,
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
    PRETRAINED_TOKENIZER_NAME_CLS_MAP,
)

check_min_version("4.39.3")
logger = logging.getLogger(__name__)

def main():
    parser = HfArgumentParser((
        OmicRawDataConfig, 
        OmicPretrainedModelAndTokenizationConfig, 
        OmicFormerTrainingArguments))
    raw_ds_config, pretrained_model_tokenizer_config, training_config = parser.parse_args_into_dataclasses()

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
    # load dataset
    #################################################################
    # load seq dataset
    seq_hf_ds = load_from_disk(raw_ds_config.raw_seq_dataset_dir)
    logger.info(f"seq HF datasets: \n{seq_hf_ds}")
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

    # scrna_raw_ds = sc.read_h5ad(raw_ds_config.raw_scrna_file_name)

    #################################################################
    # load tokenizer
    #################################################################
    # >>> for SEQ data

    if pretrained_model_tokenizer_config.seq_model_name == PretrainedModelName.HYENADNA:
        tokenizer_cls = PRETRAINED_TOKENIZER_NAME_CLS_MAP[pretrained_model_tokenizer_config.seq_model_name]
        tokenizer = tokenizer_cls(
            characters=['A', 'C', 'G', 'T', 'N'],  # add DNA characters, 7, 8, 9, 10, 11
            model_max_length=pretrained_model_tokenizer_config.seq_model_input_max_len
        )
    else:
        raise ValueError(f"Unsupported pretrained model name: {pretrained_model_tokenizer_config.seq_model_name}")

    logger.info(f"Sequence pretrained model {pretrained_model_tokenizer_config.seq_model_name.value} from "
                + f"{pretrained_model_tokenizer_config.seq_model_path}, \n  tokenizer: \n{tokenizer}")
    logger.info(f"vocab size: {tokenizer.vocab_size}")

    def seq_tokenize_fn(untokenized_ds):
        token_enocdes = tokenizer(
            untokenized_ds["seq"],
            truncation=False, 
            padding=False, # we did not padding during tokenization
            add_special_tokens=False,
            return_attention_mask=True,
        )
        return {
            "input_ids": token_enocdes["input_ids"],
            # "token_type_ids": token_type_batch,
            "attention_mask": token_enocdes["attention_mask"],
        }
    
    # NOTE: filter out the sequences with length < pretrained_model_tokenizer_config.seq_min_len
    def filter_fn(tokenized_dataset):
        return_vals = []
        for attention_mask in tokenized_dataset["attention_mask"]:
            return_vals.append(
                    sum(attention_mask) >= pretrained_model_tokenizer_config.seq_min_len
            )
        return return_vals
    
    with training_config.main_process_first(desc="SEQ dataset map tokenization"):
        _tokenized_ds = seq_hf_ds.map(
            seq_tokenize_fn,
            batched=True, 
            load_from_cache_file=True if not pretrained_model_tokenizer_config.use_streaming else None, 
            num_proc=15 if not pretrained_model_tokenizer_config.use_streaming else None,
            desc="Running tokenizer on raw SEQ dataset (map), " + "no streaming" if not pretrained_model_tokenizer_config.use_streaming else "use streaming",
        ).filter(
            filter_fn,
            batched=True,
            load_from_cache_file=True if not pretrained_model_tokenizer_config.use_streaming else None,
            num_proc=10 if not pretrained_model_tokenizer_config.use_streaming else None,
            desc="Running tokenizer on raw SEQ dataset (filter), " + "no streaming" if not pretrained_model_tokenizer_config.use_streaming else "use streaming",
        )

        if pretrained_model_tokenizer_config.use_streaming:
            # this will take a while. (~6 min per chrosome)
            dd = DatasetDict()
            for ds_name, iterable_ds in _tokenized_ds.items():
                ds = Dataset.from_generator(lambda: (yield from iterable_ds), features=iterable_ds.features)
                dd[ds_name] = ds
            dd.save_to_disk(pretrained_model_tokenizer_config.tokenized_seq_dataset_dir)
            dd.cleanup_cache_files()
        else:
            _tokenized_ds.save_to_disk(pretrained_model_tokenizer_config.tokenized_seq_dataset_dir + f"/{raw_ds_config.processed_cell_type_name}")
            _tokenized_ds.cleanup_cache_files()
            # >>> SEQ tokenized data:
            # DatasetDict({
            #     train: Dataset({
            #         features: ['sample', 'pos', 'seq', 'peak_value', 'input_ids', 'attention_mask'],
            #         num_rows: 124167712
            #     })
            #     test: Dataset({
            #         features: ['sample', 'pos', 'seq', 'peak_value', 'input_ids', 'attention_mask'],
            #         num_rows: 4208429
            #     })
            # })

if __name__=='__main__':
    main()


    



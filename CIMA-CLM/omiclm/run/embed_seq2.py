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
import os
import sys
from typing import Any, Dict, List, Mapping, NewType
import logging
from pathlib import Path
import json
import gzip

import numpy as np
import torch
import torch.nn.functional as F
import transformers
import datasets
from tqdm import tqdm
from dataclasses import dataclass
 

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

InputDataClass = NewType("InputDataClass", Any)

class SeqTokenizedDataCollator:
    def __init__(self, pad_token_id = 4):
        self.pad_token_id = pad_token_id
    
    # adapted from `transformers.DefaultDataCollator`
    def __call__(self, features: List[InputDataClass]):
        if not isinstance(features[0], Mapping):
            features = [vars(f) for f in features]
        first = features[0]
        batch = {}

        # get the max length in the batch
        max_len = max([len(f['input_ids']) for f in features]) 

        # Handling of all other possible keys.
        # Again, we will use the first element to figure out which key/values are not None for this model.
        for k, v in first.items():
            if v is not None and not isinstance(v, str):
                batch[k] = torch.tensor([f[k] + (max_len - len(f[k])) * [self.pad_token_id] 
                                         if k=="input_ids" else f[k] for f in features])
        return batch


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
    # load TOKENIZED dataset
    #################################################################

    tokenized_dataset_dir = pretrained_model_tokenizer_config.tokenized_seq_dataset_dir
    logger.info(f"loading seq tokenized dataset from: {tokenized_dataset_dir}")

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

        train_val_dataset = seq_hf_ds['train'].train_test_split(test_size=0.1, shuffle=True, seed=42)
        
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

    # #################################################################
    # # load model
    # #################################################################
    seq_model_dir = pretrained_model_tokenizer_config.seq_model_path
    logger.info(f"loading sequence pretrained model {pretrained_model_tokenizer_config.seq_model_name.value} from \n"
                + f"{pretrained_model_tokenizer_config.seq_model_path}.")
    seq_model = PRETRAINED_MODEL_NAME_CLS_MAP[pretrained_model_tokenizer_config.seq_model_name].from_pretrained(
        seq_model_dir,
        download=False,
        config=None,
        use_head=False,
    )
    seq_model.eval()

    seq_data_collator = SeqTokenizedDataCollator(pad_token_id=4)
    t_args = TrainingArguments(
        output_dir = training_config.output_dir + '/seq_embedding',
        do_train = False,
        do_predict = True,
        per_device_eval_batch_size = 192,    
        dataloader_drop_last = False,
        report_to=None,
        dataloader_num_workers = 8,
        disable_tqdm=True,
        eval_accumulation_steps=16, # see https://discuss.huggingface.co/t/batch-size-for-trainer-predict/3374
        ddp_timeout=60 * 60 * 2,
    )
    trainer = Trainer(model = seq_model, args = t_args, data_collator=seq_data_collator)

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
    
    n_shard = 5000 # 27937735
    for _dstype, _token_ds in zip(dstype, tokenized_ds):
        emb_save_dir = os.path.join(pretrained_model_tokenizer_config.emb_seq_dataset_dir, _dstype)
        if not os.path.exists(emb_save_dir):
            os.makedirs(emb_save_dir, exist_ok=True)

        for idx_shard in range(n_shard):

            with training_config.main_process_first(desc="chunking data"):
                i_dataset = _token_ds.shard(num_shards=n_shard, index=idx_shard)
                i_dataset = i_dataset.add_column('row_index', list(range(len(i_dataset))))

            logger.info(f">>>embedding {idx_shard}/{n_shard} shard with {len(i_dataset)} samples..")

            _pred = trainer.predict(i_dataset)

            record_ids = _pred.predictions[0]

            assert len(record_ids) == len(i_dataset), f'predicted number of records n={len(record_ids)} not equal to the number of sharded i_dataset n={len(i_dataset)}'

            logger.info(f">>>post-process and save predicted reults {idx_shard}/{n_shard}.")
            with training_config.main_process_first(desc="post-processing and save data"):
                emb_list = []
                for i in range(len(record_ids)):
                    _record_id = int(record_ids[i])
                    _row_index = int(_pred.predictions[1][i])
                    example = i_dataset[_row_index]
                    sample = example['sample']
                    pos = example['pos']
                    peak_value = example['peak_value']
                    _len = len(example['seq'])
                    
                    # seq_embedding = _pred.predictions[2][i][:_len,:].tolist() # result file too large in this way.
                    
                    # adaptive avg pooling
                    emb = F.adaptive_avg_pool1d(torch.tensor(_pred.predictions[2][i][:_len,:]).T, 20).T
                    emb = emb.numpy().tolist()

                    # select correct data...
                    emb_list.append(
                        {
                            'record_id': _record_id,
                            'sample':sample,
                            'pos': pos,
                            'seq_len': _len,
                            'peak_value': peak_value,
                            'seq_embedding':emb,
                        }
                    )

                # i_file_name = os.path.join(emb_save_dir, f"{_dstype}_{idx_shard}.json.gz")
                # logger.info(f">>>saving {idx_shard}/{n_shard} chunks at {i_file_name}")
                # with gzip.open(i_file_name, 'wt', encoding="utf-8") as emb_f:
                #     json.dump(emb_list, emb_f)
            with training_config.main_process_first(desc="post-processing and save data"):      
                Dataset.from_list(emb_list).save_to_disk(
                    os.path.join(emb_save_dir, f"{_dstype}_{idx_shard}"),
                    num_proc=64
                )

                # for debug
                # if idx_shard == 3:
                #     break

        logger.info(f"clean up {_dstype} cache files.")
        with training_config.main_process_first(desc="cleaning up cache files"):
            _token_ds.cleanup_cache_files()
    
    logger.info("clean up train_test_split cache files.")
    with training_config.main_process_first(desc="cleaning up cache files"):
        train_val_dataset.cleanup_cache_files()
    logger.info("DONE")


if __name__=='__main__':
    main()


    



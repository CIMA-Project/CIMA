#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
@File    :		_trainer.py
@Time    :   	2024/06/10 23:10:08
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
import torch
import torch.nn as nn
from transformers import Trainer
import time

from ._loss import MSLELoss, MSEMSLELoss


class OmicFormerTrainer(Trainer):
    
    def __init__(
        self,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
    
    def compute_loss(self, model, inputs, return_outputs=False):
        """
        How the loss is computed by Trainer. By default, all models return the loss in the first element.
        """
        # t_start = time.time()
        peak_value = inputs.pop("peak_value")
        # loss_fct = MSEMSLELoss(alpha=self.args.msle_alpha, beta=self.args.msle_beta, auxiliary_loss=self.args.auxiliary_loss)
        # NOTE we use MSELoss for log1p_norm_peak_value
        if self.args.loss_fn_name == "mse":
            loss_fct = nn.MSELoss()
        elif self.args.loss_fn_name == "msle_mse":
             loss_fct = MSEMSLELoss(alpha=self.args.msle_alpha, beta=self.args.msle_beta, auxiliary_loss=self.args.auxiliary_loss)
        elif self.args.loss_fn_name == "msle":
            loss_fct = MSLELoss() 
        else:
            raise ValueError(f"loss_fn_name: {self.args.loss_fn_name} is not supported! Using 'mse', 'msle_mse' or 'msle'.")
        outputs = model(**inputs)
        loss = loss_fct(outputs['logits'].view(-1, self.model.config.n_outputs).contiguous(), 
                        peak_value.view(-1, self.model.config.n_outputs).contiguous())
        
        # torch.cuda.current_stream().synchronize()  
        # t_seq_emd_end = time.time()
        # print(f"LOSS time: {t_seq_emd_end - t_start}")

        # if running on single GPU with DDP (multi GPU setup), then the following error will be raised:
        #   RuntimeError: Expected to have finished reduction in the prior iteration before starting a new one. 
        #   This error indicates that your module has parameters that were not used in producing loss. 
        #   You can enable unused parameter detection by passing the keyword argument `find_unused_parameters=True` to `torch.nn.parallel.DistributedDataParallel`, and by 
        #   making sure all `forward` function outputs participate in calculating loss. 
        #   If you already have done that, then the distributed data parallel module wasn't able to locate the output 
        #   tensors in the return value of your module's `forward` function. Please include the loss function and the structure 
        #   of the return value of `forward` of your module when reporting this issue (e.g. list, dict, iterable).
        #   See: https://github.com/pytorch/pytorch/issues/43259#issuecomment-706486925
        return (loss, outputs) if return_outputs else loss
    
    def save_model(self, output_dir, _internal_call:bool = False):
        self.model.save_pretrained(output_dir, safe_serialization=False)
        # self.tokenizer.save_pretrained(output_dir)
    


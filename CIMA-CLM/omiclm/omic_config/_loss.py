#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
@File    :		_loss.py
@Time    :   	2024/06/07 09:48:45
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
import torch.nn.functional as F


class MSLELoss(nn.Module):
    """
    Mean Squared Logarithmic Error

    """
    def __init__(self) -> None:
        super().__init__()
        self.mse = nn.MSELoss()
    
    def forward(self, input: torch.Tensor, target: torch.Tensor) -> torch.Tensor:

        msle = self.mse(torch.log1p(input), torch.log1p(target))

        return msle
    

class MSEMSLELoss(nn.Module):
    """
    Mean Squared Logarithmic Error

    """
    def __init__(self, alpha=2.0, beta=10.0, auxiliary_loss='mse') -> None:
        super().__init__()
        self.alpha = alpha
        self.beta = beta

        self.msle_loss = MSLELoss()

        if auxiliary_loss == 'mse':
            self.au_loss = nn.MSELoss()
        else:
            self.au_loss = nn.L1Loss()
        
    
    def forward(self, input: torch.Tensor, target: torch.Tensor) -> torch.Tensor:

        loss = self.alpha * self.au_loss(input, target) + self.beta * self.msle_loss(input, target)
        # loss = self.alpha * au_loss
        return loss

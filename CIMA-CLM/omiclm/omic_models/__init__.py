#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
@File    :		__init__.py
@Time    :   	2024/05/30 11:56:42
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

from .pretrained_model import (
    HyenaDNAPreTrainedModel,
    GeneformerModelWrapper,
    ScGPTModelWrapper,

    HyenaDNACharacterTokenizer,
    GeneformerTokenizerWrapper,
    ScGPTTokenizerWrapper,

    GeneformerDataProcessor,
    ScGPTDataProcessor,

    ExplicitEnum,
    PretrainedModelName,

    PRETRAINED_MODEL_NAME_CLS_MAP,
    PRETRAINED_TOKENIZER_NAME_CLS_MAP,
    SCRNA_DATA_PREPROCESSOR_FOR_PRETRAINED_MODEL_MAP

)

from .omic_model import (
    OmicFormerConfig,
    OmicFormerPreTrainedModel,
    OmicFormer,
)
# from dataclasses import dataclass
# from enum import Enum
# from typing import Optional

# from .pretrained_model import (
#     HyenaDNAPreTrainedModel,
#     GeneformerModelWrapper,
#     ScGPTModelWrapper,

#     HyenaDNACharacterTokenizer,
#     GeneformerTokenizerWrapper,
#     ScGPTTokenizerWrapper,

#     GeneformerDataProcessor,
#     ScGPTDataProcessor
# )

# class ExplicitEnum(Enum):
#     """
#     Enum with more explicit error message for missing values.
#     """

#     @classmethod
#     def _missing_(cls, value):
#         raise ValueError(
#             f"{value} is not a valid {cls.__name__}, please select one of {list(cls._value2member_map_.keys())}"
#         )

# class PretrainedModelName(ExplicitEnum):
#     HYENADNA = "hyenadna"
#     GENEFORMER = "geneformer"
#     SCGPT = "scgpt"
#     SCVI = "scvi"


# @dataclass
# class OmicLMModelInput:
#     sample_id: Optional[str] = None,
#     cell_type: Optional[str] = None,
#     input_ids=None


# PRETRAINED_MODEL_NAME_CLS_MAP = {
#     PretrainedModelName.HYENADNA: HyenaDNAPreTrainedModel,
#     PretrainedModelName.GENEFORMER: GeneformerModelWrapper,
#     PretrainedModelName.SCGPT: ScGPTModelWrapper,
# }

# PRETRAINED_TOKENIZER_NAME_CLS_MAP = {
#     PretrainedModelName.HYENADNA: HyenaDNACharacterTokenizer,
#     PretrainedModelName.GENEFORMER: GeneformerTokenizerWrapper,
#     PretrainedModelName.SCGPT: ScGPTTokenizerWrapper,
# }

# SCRNA_DATA_PREPROCESSOR_FOR_PRETRAINED_MODEL_MAP = {
#     PretrainedModelName.GENEFORMER: GeneformerDataProcessor,
#     PretrainedModelName.SCGPT: ScGPTDataProcessor
# }




from dataclasses import dataclass
from enum import Enum
from typing import Optional

from .hyenadna import HyenaDNAPreTrainedModel, HyenaDNAModel, HyenaDNACharacterTokenizer

from .scgpt import ScGPTConfig, ScGPTModelWrapper, ScGPTTokenizerWrapper

from .geneformer import GeneformerModelWrapper, GeneformerTokenizerWrapper

from .data_prep import GeneformerDataProcessor, ScGPTDataProcessor


class ExplicitEnum(Enum):
    """
    Enum with more explicit error message for missing values.
    """

    @classmethod
    def _missing_(cls, value):
        raise ValueError(
            f"{value} is not a valid {cls.__name__}, please select one of {list(cls._value2member_map_.keys())}"
        )

class PretrainedModelName(ExplicitEnum):
    HYENADNA = "hyenadna"
    GENEFORMER = "geneformer"
    SCGPT = "scgpt"
    SCVI = "scvi"

PRETRAINED_MODEL_NAME_CLS_MAP = {
    PretrainedModelName.HYENADNA: HyenaDNAPreTrainedModel,
    PretrainedModelName.GENEFORMER: GeneformerModelWrapper,
    PretrainedModelName.SCGPT: ScGPTModelWrapper,
}

PRETRAINED_TOKENIZER_NAME_CLS_MAP = {
    PretrainedModelName.HYENADNA: HyenaDNACharacterTokenizer,
    PretrainedModelName.GENEFORMER: GeneformerTokenizerWrapper,
    PretrainedModelName.SCGPT: ScGPTTokenizerWrapper,
}

SCRNA_DATA_PREPROCESSOR_FOR_PRETRAINED_MODEL_MAP = {
    PretrainedModelName.GENEFORMER: GeneformerDataProcessor,
    PretrainedModelName.SCGPT: ScGPTDataProcessor
}

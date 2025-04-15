#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
@File    :		_scgpt_model.py
@Time    :   	2024/05/26 14:57:06
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
import json
import logging
from typing import Any, List, Mapping, Optional, Union
import numpy as np

import torch
from torchtext.vocab import Vocab

from anndata import AnnData
import scanpy as sc

from transformers import PretrainedConfig, PreTrainedModel
from transformers.utils import WEIGHTS_NAME, CONFIG_NAME
from transformers.utils.hub import cached_file

from scgpt.model import TransformerModel
from scgpt.tokenizer.gene_tokenizer import GeneVocab

from scgpt import tasks

logger = logging.getLogger(__name__)

class ScGPTConfig(PretrainedConfig):

    model_type = 'scgpt_transformer'

    def __init__(
        self,
        ntoken: int = 60697,
        d_model: int = 512,
        nhead: int = 8,
        d_hid: int = 512,
        nlayers: int = 12,
        nlayers_cls: int = 3,
        n_cls: int = 1,
        # vocab: Any = None,
        dropout: float = 0.2,
        pad_token: str = "<pad>",
        pad_value: int = -2,
        do_mvc: bool = True,
        do_dab: bool = False,
        use_batch_labels: bool = False,
        num_batch_labels: Optional[int] = None,
        domain_spec_batchnorm: Union[bool, str] = False,
        input_emb_style: str = "continuous",
        n_input_bins: Optional[int] = 51,
        cell_emb_style: str = "cls",                        # could be one of "cls", "avg-pool", "w-pool"
        mvc_decoder_style: str = "inner product",
        ecs_threshold: float = 0.3,
        explicit_zero_prob: bool = False,
        use_fast_transformer: bool = True,
        fast_transformer_backend: str = "flash",
        pre_norm: bool = False,
        mask_value: int = -1,
        embsize: int = 512,
        **kwargs,
    ):
        """
        Paramter are initialized with the value from args.json in the scGPT_human pretrained model.
        
        See:
        https://github.com/bowang-lab/scGPT

        """
        self.ntoken = ntoken
        self.d_model = d_model
        self.nhead = nhead 
        self.d_hid = d_hid 
        self.nlayers = nlayers 
        self.nlayers_cls = nlayers_cls                                  
        self.n_cls = n_cls
        # self.vocab = vocab                                            # assign separately                   
        self.dropout = dropout
        self.pad_token = pad_token
        self.pad_value = pad_value
        self.do_mvc = do_mvc
        self.do_dab = do_dab
        self.use_batch_labels = use_batch_labels
        self.num_batch_labels = num_batch_labels
        self.domain_spec_batchnorm = domain_spec_batchnorm
        self.input_emb_style = input_emb_style
        self.n_input_bins = n_input_bins
        self.cell_emb_style = cell_emb_style
        self.mvc_decoder_style = mvc_decoder_style
        self.ecs_threshold= ecs_threshold
        self.explicit_zero_prob = explicit_zero_prob
        self.use_fast_transformer = use_fast_transformer
        self.fast_transformer_backend = fast_transformer_backend
        self.pre_norm = pre_norm
        self.mask_value = mask_value
        self.embsize = embsize

        super().__init__(architectures=["TransformerModel"], **kwargs)


class ScGPTModelWrapper(PreTrainedModel):

    config_class = ScGPTConfig

    def __init__(
        self,
        config: ScGPTConfig,
        vocab: Optional[Vocab] = None, 
        device = None,
        **kwargs,
    ) -> None:

        super().__init__(config, **kwargs)
        
        self.vocab = vocab

        # NOTE: scGPT does not allow dtype config
        model_args = dict(
            ntoken=config.ntoken,
            d_model=config.d_model,
            nhead=config.nhead,
            d_hid=config.d_hid,
            nlayers=config.nlayers,
            nlayers_cls=config.nlayers_cls,
            n_cls=config.n_cls,
            dropout=config.dropout,
            pad_token=config.pad_token,
            pad_value=config.pad_value,
            do_mvc=config.do_mvc,
            do_dab=config.do_dab,
            use_batch_labels=config.use_batch_labels,
            domain_spec_batchnorm=config.domain_spec_batchnorm,
            explicit_zero_prob=config.explicit_zero_prob,
            use_fast_transformer=config.use_fast_transformer,
            fast_transformer_backend=config.fast_transformer_backend,
            pre_norm=config.pre_norm,
        )
        self.model = TransformerModel(vocab=vocab,**model_args)
        self.model.to(device)
    
    def forward(self, **kwargs):
        # TODO:
        msg = "TODO, current scGPT is only used for embedding." 
        logger.warning(msg)
        return self.model(**kwargs)

    # adapted from scgpt.tasks.cell_emb.emb_data
    def extract_sample_embedding(
        self,
        adata_or_file,
        gene_col: str="gene_col",
        max_length=1200,
        cell_embedding_mode:str="cls",      #    could be one of "cls", "avg-pool", "w-pool"   
        batch_size: int = 64, 
        embedding_key: str = "cell_embedding",
        obs_to_save: Optional[List[str]] = "sample",
        return_new_adata=True, 
        **kwargs,     
    ):
        """ embeddimg cells by the pretrained scGPT

        adata_or_file - h5ad data should at least be processed:
            # >>> preprocess, See `data_prep._scgpt_data_processor.py`
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            # highly variable genes
            sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor='seurat_v3')
            adata = adata[:, adata.var['highly_variable']]

        Parameters
        ----------
        adata_or_file : _type_
            _description_
        gene_col : str, optional
            _description_, by default "gene_col"
        max_length : int, optional
            _description_, by default 1200
        cell_embedding_mode : str, optional
            _description_, by default "cls"
        obs_to_save : Optional[List[str]], optional
            _description_, by default None
        return_new_adata : bool, optional
            _description_, by default True

        Returns
        -------
        _type_
            _description_
        """
        
        if isinstance(adata_or_file, AnnData):
            adata = adata_or_file
        else:
            adata = sc.read_h5ad(adata_or_file)
        if isinstance(obs_to_save, str):
            assert obs_to_save in adata.obs, f"obs_to_save {obs_to_save} not in adata.obs"
            obs_to_save = [obs_to_save]

        # verify gene col
        if gene_col == "index":
            adata.var["index"] = adata.var.index
        else:
            assert gene_col in adata.var
        
        pad_token = "<pad>"
        special_tokens = [pad_token, "<cls>", "<eoc>"]
        for s in special_tokens:
            if s not in self.vocab:
                self.vocab.append_token(s)

        adata.var["id_in_vocab"] = [
            self.vocab[gene] if gene in self.vocab else -1 for gene in adata.var[gene_col]
        ]

        gene_ids_in_vocab = np.array(adata.var["id_in_vocab"])
        logger.info(
            f"match {np.sum(gene_ids_in_vocab >= 0)}/{len(gene_ids_in_vocab)} genes "
            f"in vocabulary of size {len(self.vocab)}."
        )

        adata = adata[:, adata.var["id_in_vocab"] >= 0]

        genes = adata.var[gene_col].tolist()
        gene_ids = np.array(self.vocab(genes), dtype=int)

        self.model.eval()

        # Binning will be applied after tokenization. A possible way to do is to use the unified way of binning in the data collator.
        # get cell embeddings
        # NOTE: only support cell_embedding_mode="cls"
        # For others, see TransformerModel._get_cell_emb_from_layer
        cell_embeddings = tasks.get_batch_cell_embeddings(
            adata,
            cell_embedding_mode=cell_embedding_mode,
            model=self.model,
            vocab=self.vocab,
            max_length=max_length,
            batch_size=batch_size,
            model_configs=self.config.to_dict(),
            gene_ids=gene_ids,
            use_batch_labels=False,
        )

        if return_new_adata:
            obs_df = adata.obs[obs_to_save] if obs_to_save is not None else None
            return sc.AnnData(X=cell_embeddings, obs=obs_df, dtype="float32")

        adata.obsm[embedding_key] = cell_embeddings
        return adata



    @classmethod
    def from_pretrained(
        cls,
        pretrained_model_name_or_path: str,
        config_file_name: str="args.json",
        vocab_file_name: str="vocab.json", 
        model_file_name: str="best_model.pt",
        local_files_only: bool = True,
        device=None,
        **kwargs,
    ):
        config_fname = os.path.join(pretrained_model_name_or_path, config_file_name)
        vocab_fname = os.path.join(pretrained_model_name_or_path, vocab_file_name)
        
        if not os.path.isfile(vocab_fname):
            raise FileNotFoundError(f"vocab file not found: {vocab_fname}")
        
        # read vocab
        scgpt_vocab = GeneVocab.from_file(vocab_fname)
        scgpt_vocab.set_default_index(scgpt_vocab["<pad>"])

        if not os.path.exists(config_fname):
            logger.warning(f"config file {config_fname} not found. Using default model config instead")
            config = ScGPTConfig(ntoken=len(scgpt_vocab))
        elif config_file_name == "args.json":
            model_configs = load_config_hf(
                pretrained_model_name_or_path, 
                file_name=config_file_name,
                local_files_only=local_files_only)
            try:
                d_model=model_configs["embsize"]
            except:
                d_model=model_configs["d_model"]
            try: 
                nlayers_cls = model_configs["n_layers_cls"]
            except:
                nlayers_cls = model_configs["nlayers_cls"]
            try:
                use_fast_transformer=model_configs["fast_transformer"]
            except:
                use_fast_transformer=model_configs["use_fast_transformer"]

            config = ScGPTConfig(
                ntoken=len(scgpt_vocab),
                d_model=d_model,
                nhead=model_configs["nheads"],
                d_hid=model_configs["d_hid"],
                nlayers=model_configs["nlayers"],
                nlayers_cls=nlayers_cls,
                n_cls=1,
                dropout=model_configs["dropout"],
                pad_token=model_configs["pad_token"],
                pad_value=model_configs["pad_value"],
                do_mvc=True,
                do_dab=False,
                use_batch_labels=False,
                domain_spec_batchnorm=False,
                explicit_zero_prob=False,
                use_fast_transformer=use_fast_transformer,
                fast_transformer_backend="flash",
                pre_norm=False,
                mask_value=model_configs["mask_value"],
                embsize = d_model,
            )
        elif config_file_name == "scgpt_config.json":
            config_dict = load_config_hf(
                pretrained_model_name_or_path, 
                file_name=config_file_name,
                local_files_only=local_files_only)
            config = ScGPTConfig(**config_dict)

        model = cls(config, vocab=scgpt_vocab, device=device, **kwargs)
        model_fname = os.path.join(pretrained_model_name_or_path, model_file_name)
        load_pretrained(
            model.model, 
            torch.load(model_fname),
            verbose=False
        )
        # model.model.load_state_dict(
        #     load_state_dict_hf(
        #         pretrained_model_name_path, 
        #         file_name=model_file_name, 
        #         local_files_only=local_files_only, 
        #         device=device
        #     )
        # )

        return model
    
    def save_pretrained(self, save_directory,
                        congif_file_name="scgpt_config.json", 
                        model_file_name="scgpt_pytorch_model.bin"):
        """
        Minimal implementation of save_pretrained for MambaLMHeadModel.
        Save the model and its configuration file to a directory.
        """
        # Ensure save_directory exists
        if not os.path.exists(save_directory):
            # multiple GPUs will raise FileExistsError if exist_ok=False(default value)
            os.makedirs(save_directory, exist_ok=True)

        # Save the model's state_dict
        model_path = os.path.join(save_directory, model_file_name)
        torch.save(self.state_dict(), model_path)

        # Save the configuration of the model
        config_path = os.path.join(save_directory, congif_file_name)
        self.config.to_json_file(config_path)


# copy from scgpt.utils.util
def load_pretrained(
    model: torch.nn.Module,
    pretrained_params: Mapping[str, torch.Tensor],
    strict: bool = False,
    prefix: Optional[List[str]] = None,
    verbose: bool = True,
) -> torch.nn.Module:
    """
    Load pretrained weights to the model.

    Args:
        model (torch.nn.Module): The model to load weights to.
        pretrained_params (Mapping[str, torch.Tensor]): The pretrained parameters.
        strict (bool): Whether to strictly enforce that the keys in :attr:`pretrained_params`
            match the keys returned by this module's :meth:`Module.state_dict`. Default to False.
        prefix (List[str]): The list of prefix strings to match with the keys in
            :attr:`pretrained_params`. The matched keys will be loaded. Default to None.

    Returns:
        torch.nn.Module: The model with pretrained weights.
    """

    use_flash_attn = getattr(model, "use_fast_transformer", True)
    if not use_flash_attn:
        pretrained_params = {
            k.replace("Wqkv.", "in_proj_"): v for k, v in pretrained_params.items()
        }

    if prefix is not None and len(prefix) > 0:
        if isinstance(prefix, str):
            prefix = [prefix]
        pretrained_params = {
            k: v
            for k, v in pretrained_params.items()
            if any(k.startswith(p) for p in prefix)
        }

    model_dict = model.state_dict()
    if strict:
        if verbose:
            for k, v in pretrained_params.items():
                logger.info(f"Loading parameter {k} with shape {v.shape}")
        model_dict.update(pretrained_params)
        model.load_state_dict(model_dict)
    else:
        if verbose:
            for k, v in pretrained_params.items():
                if k in model_dict and v.shape == model_dict[k].shape:
                    logger.info(f"Loading parameter {k} with shape {v.shape}")
        pretrained_params = {
            k: v
            for k, v in pretrained_params.items()
            if k in model_dict and v.shape == model_dict[k].shape
        }
        model_dict.update(pretrained_params)
        model.load_state_dict(model_dict)

    return model

def verify_exist(f_name):
    if not os.path.exists(f_name):
        raise FileNotFoundError(f"{f_name} not found.")


def read_json(json_fname):

    verify_exist(json_fname)
    
    with open(json_fname, "rt", encoding="utf-8") as f:
        dt_conf = json.load(f)
    return dt_conf


# adapted from from mamba_ssm.utils.hf import load_config_hf, load_state_dict_hf
def load_config_hf(
        model_name_or_dict, 
        file_name=None, 
        local_files_only=True
    ):
    resolved_archive_file = cached_file(
        model_name_or_dict, 
        filename=file_name if file_name is not None else CONFIG_NAME, 
        local_files_only=local_files_only, 
        _raise_exceptions_for_missing_entries=False)

    return read_json(resolved_archive_file)


def load_state_dict_hf(
        model_name_or_dict, 
        file_name=None, 
        local_files_only=True, 
        device=None, dtype=None
    ):
    # If not fp32, then we don't want to load directly to the GPU
    mapped_device = "cpu" if dtype not in [torch.float32, None] else device
    resolved_archive_file = cached_file(
        model_name_or_dict, 
        filename=file_name if file_name is not None else WEIGHTS_NAME, 
        local_files_only=local_files_only, 
        _raise_exceptions_for_missing_entries=False)
    state_dict = torch.load(resolved_archive_file, map_location=mapped_device)
    # Convert dtype before moving to GPU to save memory
    if dtype is not None:
        state_dict = {k: v.to(dtype=dtype) for k, v in state_dict.items()}
    state_dict = {k: v.to(device=device) for k, v in state_dict.items()}
    return state_dict

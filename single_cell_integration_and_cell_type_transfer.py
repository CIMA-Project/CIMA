"""
**************************************
*  @Author  :   Yuhui Zheng
*  @Time    :   2024/12/19
*  @Project :   CIMA
*  @FileName:   single_cell_integration_and_cell_type_transfer.py
*  @Summary :   Using scGLUE to build an integrated model of scRNA-seq and scATAC-seq and transfer cell type annotations from scRNA-seq to scATAC-seq.
**************************************
"""

import os
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams
from itertools import chain
import seaborn as sns
import pandas as pd
import scanpy as sc
import cupy as cp
import time
import rapids_singlecell as rsc
from matplotlib.pyplot import rc_context
import warnings
warnings.filterwarnings("ignore")

scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)
sc.settings.set_figure_params(dpi=150, facecolor='white', transparent=True, dpi_save=600)


rna = ad.read_h5ad('CIMA/scRNA/CIMA_scRNA_Annotation.h5ad',backed='r+') #Celltype l4 Annotated scRNA-seq Data
rna[rna.obs.groupby("sample").sample(frac = 0.5,replace=False).index].copy(filename='CIMA/scRNA/CIMA_scRNA_Annotation_sampling_0.5.h5ad')

rna_sampling = ad.read_h5ad('CIMA/scRNA/CIMA_scRNA_Annotation_sampling_0.5.h5ad')
atac = ad.read_h5ad('CIMA/scATAC/CIMA_scATAC_merged.h5ad')

scglue.data.get_gene_annotation(
    rna_sampling, gtf="CIMA/scRNA/V3_genes.gtf",
    gtf_by="gene_name"
)
rna_sampling.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head()

guidance = scglue.genomics.rna_anchored_guidance_graph(rna_sampling, atac)
scglue.graph.check_graph(guidance, [rna_sampling, atac])
atac.var.highly_variable.value_counts()

scglue.data.lsi(atac, n_components=100, n_iter=15, use_highly_variable=True)

scglue.models.configure_dataset(
    rna_sampling, "NB", use_highly_variable=True,
    use_layer="counts", use_rep="X_pca",use_batch='sample'
)

scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_lsi",use_batch='sample'
)

guidance_hvf = guidance.subgraph(chain(
    rna_sampling.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
)).copy()

glue = scglue.models.fit_SCGLUE(
    {"rna": rna_sampling, "atac": atac}, guidance_hvf,init_kws={"h_dim":512, "random_seed":666},
    fit_kws={"directory": "glue_0.5_8192","data_batch_size":8192}
)
glue.save("glue_0.5_8192.dill")

dx = scglue.models.integration_consistency(
    glue, {"rna": rna_sampling, "atac": atac}, guidance_hvf
)
_ = sns.lineplot(x="n_meta", y="consistency", data=dx).axhline(y=0.05, c="darkred", ls="--")

rna_sampling.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)

scglue.data.transfer_labels(rna_sampling, atac, "cell_type_l4", use_rep="X_glue", n_jobs=-1)

rna_sampling.obs['domain'] = 'RNA'
atac.obs['domain'] = 'ATAC'
combined = ad.concat([rna_sampling, atac])

atac.write_h5ad('CIMA/scATAC/CIMA_scATAC_Annotation_Transfered.h5ad')
combined.write_h5ad('CIMA/GLUE/CIMA_Combined.h5ad')

cell_type_l4_colors={
    'ILC2_IL2RA':'#9d5294',
    'NK_bright_XCL1':'#672368',
    'Transitional_NK_GZMK':'#402253',
    'Terminal_NK_dim_CD160neg':'#b69ec5',
    'Mature_NK_dim_FCGR3A':'#945a7b',
    'Inflamed_NK_dim_IFIT1':'#794f7f',
    'Cycling_NK_MKI67':'#9f7aaf',
    'Transitional_B_SOX4':'#f6d661',
    'Bn_TCL1A':'#f8f199',
    'Bn_IL6':'#f0bb40',
    'Bn_IFIT3':'#cfb649',
    'Atypical_Bm_ITGAX':'#f6d1e0',
    'Unswitched_Bm_IL6':'#e68d9a',
    'Unswitched_Bm_CD1C':'#eb6b87',
    'pre-Switched_Bm_IFIT3':'#d47da0',
    'pre-Switched_Bm_JAM3':'#dc7460',
    'Switched_Bm_IGHDneg':'#c62b21',
    'Switched_Bm_IGHE':'#a22124',
    'Switched_activated_Bm_CD86':'#d71718',
    'Plasmablast_MKI67':'#f49f54',
    'Plasma_IGHA1':'#ea553d',
    'Plasma_IGHG1':'#ec8300',
    'HSPC_CD34':'#00c982',
    'MK_GP9':'#908e5c',
    'cDC_CSF2RA':'#bcbd94',
    'cDC1_BATF3':'#bdbfa8',
    'cDC2_CD1C':'#99b194',
    'pDC_IRF4':'#dedecf',
    'AS_DC':'#8d6c5e',
    'cMono_CD14':'#749248',
    'cMono_CXCL10':'#736b4a',
    'cMono_IL1B':'#10383d',
    'cMono_IFI44L':'#bccc9b',
    'intMono_GFRA2':'#c5bdb0',
    'ncMono_C1QA':'#a7a685',
    'ncMono_FCGR3A':'#d1c59e',
    'ncMono_IFIT1':'#78826c',
    'CD4_Tn_SOX4':'#1f4c97',
    'CD4_Tn_CCR7':'#abcbe7',
    'CD4_Tn_CXCR5':'#282c71',
    'CD4_Tn_LIMS1':'#146192',
    'CD4_Tcm_CXCR5':'#64a2bf',
    'CD4_Tcm_IFI44L':'#1a3975',
    'CD4_Tem_CCR7neg':'#6c8ec2',
    'CD4_Tem_CCR5':'#5e93b9',
    'CD4_Tfh-like_CXCR5':'#306ca4',
    'CD4_Th_CR1':'#043b6e',
    'CD4_Th_LMNA':'#5671af',
    'CD4_Th_CCR4':'#111b39',
    'CD4_Th_TNFRSF11A':'#2b559d',
    'CD4_Th1-like_GZMK':'#add0da',
    'CD4_Th17-like_RORC':'#1a71a8',
    'CD4_Th22-like_CCR10':'#1f2a66',
    'CD4_Tr1-like_IL10':'#166490',
    'CD4_Treg_FCRL3':'#a0b7d5',
    'CD4_Treg_FOXP3':'#85a5d1',
    'CD4_CTL_GZMH':'#2663ac',
    'CD8_Tn_CCR7':'#7c9baa',
    'CD8_Tn_SOX4':'#264751',
    'CD8_Tcm_IFI44L':'#5a9e98',
    'CD8_Tem_CCR7neg':'#408386',
    'CD8_CTL_GZMK':'#a4b7af',
    'CD8_CTL_IFI44L':'#89b29e',
    'CD8_CTL_GZMB':'#c2d7d3',
    'Cycling_T_MKI67':'#6a7c72',
    'pre-T-like_CABP4':'#62b2c6',
    'MAIT_SLC4A10':'#32a4c3',
    'gdT1_TRDV1':'#4e8ea0',
    'gdT2_IL12RB2':'#72b9bc',
    'gdT2_GZMH':'#15858f',
    'gdT2_GZMK':'#51bac8',
    'NKT_IFNG':'#188ea7',
    'NKT_NCR1':'#20a5b7'
}

cell_type_l3={
    'CD4_Tn_SOX4' : 'CD4_Tn',
    'CD4_Tn_CCR7' : 'CD4_Tn',
    'CD4_Tn_CXCR5' : 'CD4_Tn',
    'CD4_Tn_LIMS1' : 'CD4_Tn',
    'CD4_Tcm_CXCR5' : 'CD4_Tcm',
    'CD4_Tcm_IFI44L' : 'CD4_Tcm',
    'CD4_Tem_CCR7neg' : 'CD4_Tem',
    'CD4_Tem_CCR5' : 'CD4_Tem',
    'CD4_Tfh-like_CXCR5' : 'CD4_helper',
    'CD4_Th_CR1' : 'CD4_helper',
    'CD4_Th_LMNA' : 'CD4_helper',
    'CD4_Th_CCR4' : 'CD4_helper',
    'CD4_Th_TNFRSF11A' : 'CD4_helper',
    'CD4_Th1-like_GZMK' : 'CD4_helper',
    'CD4_Th17-like_RORC' : 'CD4_helper',
    'CD4_Th22-like_CCR10' : 'CD4_helper',
    'CD4_Tr1-like_IL10' : 'CD4_Tr1',
    'CD4_Treg_FCRL3' : 'CD4_Treg',
    'CD4_Treg_FOXP3' : 'CD4_Treg',
    'CD4_CTL_GZMH' : 'CD4_CTL',
    'CD8_Tn_CCR7' : 'CD8_Tn',
    'CD8_Tn_SOX4' : 'CD8_Tn',
    'CD8_Tcm_IFI44L' : 'CD8_Tcm',
    'CD8_Tem_CCR7neg' : 'CD8_Tem',
    'CD8_CTL_GZMK' : 'CD8_CTL',
    'CD8_CTL_IFI44L' : 'CD8_CTL',
    'CD8_CTL_GZMB' : 'CD8_CTL',
    'Cycling_T_MKI67' : 'Cycling_T',
    'pre-T-like_CABP4' : 'pre-T-like',
    'gdT1_TRDV1' : 'gdT1',
    'gdT2_GZMK' : 'gdT2',
    'gdT2_IL12RB2' : 'gdT2',
    'gdT2_GZMH' : 'gdT2',
    'MAIT_SLC4A10' : 'MAIT',
    'NKT_NCR1' : 'NKT',
    'NKT_IFNG' : 'NKT',
    'Transitional_B_SOX4' : 'Transitional_B',
    'Bn_TCL1A' : 'Bn',
    'Bn_IL6' : 'Bn',
    'Bn_IFIT3' : 'Bn',
    'Atypical_Bm_ITGAX' : 'Atypical_Bm',
    'Unswitched_Bm_IL6' : 'Unswitched_Bm',
    'Unswitched_Bm_CD1C' : 'Unswitched_Bm',
    'pre-Switched_Bm_IFIT3' : 'pre-Switched_Bm',
    'pre-Switched_Bm_JAM3' : 'pre-Switched_Bm',
    'Switched_Bm_IGHDneg' : 'Switched_Bm',
    'Switched_Bm_IGHE' : 'Switched_Bm',
    'Switched_activated_Bm_CD86' : 'Switched_Bm',
    'Plasmablast_MKI67' : 'Plasmablast',
    'Plasma_IGHA1' : 'Plasma',
    'Plasma_IGHG1' : 'Plasma',
    'HSPC_CD34' : 'HSPC',
    'MK_GP9' : 'MK',
    'cDC_CSF2RA' : 'cDC',
    'cDC1_BATF3' : 'cDC',
    'cDC2_CD1C' : 'cDC',
    'AS_DC' : 'pDC',
    'pDC_IRF4' : 'pDC',
    'cMono_CD14' : 'cMono',
    'cMono_CXCL10' : 'cMono',
    'cMono_IL1B' : 'cMono',
    'cMono_IFI44L' : 'cMono',
    'intMono_GFRA2' : 'intMono',
    'ncMono_C1QA' : 'ncMono',
    'ncMono_FCGR3A' : 'ncMono',
    'ncMono_IFIT1' : 'ncMono',
    'ILC2_IL2RA' : 'ILC2',
    'NK_bright_XCL1' : 'CD56_bright_NK',
    'Transitional_NK_GZMK' : 'Transitional_NK',
    'Mature_NK_dim_FCGR3A' : 'Mature_NK',
    'Inflamed_NK_dim_IFIT1' : 'Mature_NK',
    'Terminal_NK_dim_CD160neg' : 'Terminal_NK',
    'Cycling_NK_MKI67' : 'Cycling_NK'
}

cell_type_l3_color =   {'CD4_Tn': '#88afd3',
    'CD8_Tn':'#004461' ,
    'CD4_T_Helper': '#648eca',
    'cMono': '#6f8451',
    'CD8_CTL':'#9ed5cd' ,
    'Mature_NK': '#945a7b' ,
    'NKT': '#0099b5',
    'gdT2': '#009d9d',
    'CD8_Tem': '#408386',
    'MAIT':'#32a4c3' ,
    'Naive_B': '#f8f081',
    'CD4_Tem': '#006e8e',
    'ncMono': '#768d80',
    'Switched_Bm': '#e22912',
    'CD4_Tcm':'#192d5e' ,
    'CD4_Treg': '#2663ac',
    'CD4_CTL': '#2663ac',
    'Terminal_NK': '#b69ec5' ,
    'Transitional_B': '#f6d661',
    'Transitional_NK': '#402253',
    'CD8_Tcm': '#5a9e98' ,
    'cDC': '#a3b494',
    'Unswitched_Bm': '#e77694',
    'pre-Switched_Bm': '#d26077',
    'pDC': '#8e8d83',
    'Atypical_Bm':'#f6d1e0' ,
    'CD56_bright_NK': '#672368',
    'Plasma': '#ee7943',
    'intMono': '#c5bdb0',
    'MK': '#8d8c5b',
    'Cycling_T': '#6a7c72',
    'gdT1': '#4e8ea0',
    'CD4_Tr1':'#166490' ,
    'pre-T-like': '#62b2c6',
    'Cycling_NK': '#9f7aa9',
    'HSPC':'#3ab37b' ,
    'Plasmablast': '#f49f54',
    'ILC2': '#9d5294' }

atac.obs['cell_type_l3']=atac.obs['cell_type_l4'].replace(cell_type_l3)

# sampling for plot
atac_sampling = sc.pp.subsample(atac, fraction=0.3, copy=True) 

import rmm
from rmm.allocators.cupy import rmm_cupy_allocator

rmm.reinitialize(
    managed_memory=True,  # Allows oversubscription
    pool_allocator=False,  # default is False
    devices=0,  # GPU device IDs to register. By default registers only GPU 0.
)
cp.cuda.set_allocator(rmm_cupy_allocator)
rsc.get.anndata_to_GPU(atac_sampling)

sc.pp.neighbors(atac_sampling, n_neighbors=30, use_rep="X_glue", random_state=300)
rsc.tl.umap(atac_sampling, min_dist=0.5, random_state=111, spread=1)

with rc_context({'figure.figsize': (8, 8)}):
    sc.pl.umap(atac_sampling, color=['cell_type_l4'],palette = cell_type_l4_colors,
               legend_fontsize=8, legend_fontoutline=2, sort_order=True, size=1, edges_color=None, save='_CIMA_scATAC_celltype_l4.pdf')
    
with rc_context({'figure.figsize': (8, 8)}):
    sc.pl.umap(atac_sampling, color=['cell_type_l3'],palette = cell_type_l3_color,
               legend_fontsize=8, legend_fontoutline=2, sort_order=True, size=1, edges_color=None, save='_CIMA_scATAC_celltype_l3.pdf')

# sampling for plot
combined_sampling = sc.pp.subsample(combined,fraction=0.1,copy=True) 
rsc.get.anndata_to_GPU(combined_sampling)
sc.pp.neighbors(combined_sampling, n_neighbors=30, use_rep="X_glue",random_state=300)
rsc.tl.umap(combined_sampling, min_dist=0.6, random_state=5, spread=1)

sc.pl.umap(combined_sampling, color=["domain"], palette={'ATAC':'#1C3030','RNA':'#7F9B54'})

atac_sampling.write_h5ad('CIMA/scATAC/CIMA_scATAC_Annotation_Transfered_sampling.h5ad')
combined_sampling.write_h5ad('CIMA/GLUE/CIMA_Combined_sampling.h5ad')
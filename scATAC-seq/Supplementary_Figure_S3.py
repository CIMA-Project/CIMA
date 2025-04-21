#!/usr/bin/env python
# coding: utf-8

# In[1]:
import numpy as np
import pandas as pd
import anndata as ad
import os
import seaborn as sns
#
import scanpy as sc
sc.settings.set_figure_params(dpi=80, figsize=(4, 4),facecolor='white')
sns.set_style(style='ticks')
#
import matplotlib.pyplot as plt
get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'retina'")
# * loading scATAC
# In[2]:
gene_matrix=sc.read_h5ad('../scATAC/SnapATAC2/final_filtered/pbmc_final_genescore_GLUE.h5ad')
#
marker_genes = {
    "list1": ['CD3D', 'CD4', 'CD8A', 'CD79A', 'CD68'],
    "list2": ["MS4A1","NEIL1","TCL1A","CD27","XBP1"],
    "list3": ["KIR2DL1", "NCR1", "IL1R1", 'KLRD1'],
    "list4": ["CCR7", "GPR183", "CTLA4", "GZMH"],
    "list5": ["CD3D", "CD8B", "CCR7", "GPR183",'GZMH'],
    "list6": ["MKI67", "NCR1", "SLC4A10", "TRDC"],
    "list7": ["CD34", "FCN1", "HLA-DPB1", "LILRA4",'PPBP'],
    "list8": ["IL7R", "KIT", "NKG7", "NCAM1", 'B3GAT1', 'FCGR3A', 'STMN1'],
}
marker_genes_in_data = {key: [i for i in marker_genes[key] if i in gene_matrix.var_names.values]
                for key in marker_genes.keys()
               }
# In[3]:
gene_matrix
# * visualization
# In[5]:
sc.settings.set_figure_params(dpi=80, figsize=(4, 4),facecolor='white')
#
for ct in marker_genes_in_data.keys():
    print(f"{ct.upper()}:")  # print cell subtype name
    sc.pl.umap(
        gene_matrix,
        color=marker_genes_in_data[ct],
        vmin=0,
        vmax=4,  # set vmax to the 99th percentile of the gene count instead of the maximum, to prevent outliers from making expression in other cells invisible. Note that this can cause problems for extremely lowly expressed genes.
        sort_order=False,  # do not plot highest expression on top, to not get a biased view of the mean expression among cells
        frameon=False,
        cmap="YlOrRd",  # or choose another color map e.g. from here: https://matplotlib.org/stable/tutorials/colors/colormaps.html
        #save=f'_scATAC_markers_{ct}.pdf'
    )
    print("\n\n\n")  # print white space for legibility


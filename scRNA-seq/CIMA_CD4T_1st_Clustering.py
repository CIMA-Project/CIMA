import os
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
import scanpy as sc
import anndata
import cupy as cp
import numpy as np
import pandas as pd
import time
import rapids_singlecell as rsc
from matplotlib.pyplot import rc_context
from rapids_singlecell.cunnData import cunnData

import warnings
warnings.filterwarnings("ignore")

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=150, facecolor='white') 

import rmm
rmm.reinitialize(
    managed_memory=True)
cp.cuda.set_allocator(rmm.rmm_cupy_allocator)

adata = sc.read_h5ad('CIMA/scRNA/CIMA_CD4T_doublet_2.4M.h5ad')
adata.X = adata.layers['counts'].copy()

sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata,batch_key='sample')
adata.var.highly_variable.value_counts()

mt_genes = list(adata.var_names[adata.var_names.str.match(r'^MT-')])
rp_genes = list(adata.var_names[adata.var_names.str.match(r'^RP[SL]')]) 
ncRNA_genes = list(adata.var_names[adata.var_names.str.match(r'^[A-Z][A-Z][0-9].*\.[0-9]')])
LINC_genes = list(adata.var_names[adata.var_names.str.match(r'(^LOC|LINC)[1-9]*')])
remove_genes = mt_genes + rp_genes + ncRNA_genes + LINC_genes

hvg_list = adata.var_names[adata.var.highly_variable].to_list()
pd.DataFrame(hvg_list).to_excel('./CD4T_1st_hvg_list.xlsx')
hvg_remove = [gene for gene in hvg_list if gene in remove_genes]
pd.DataFrame(hvg_remove).to_excel('./CD4T_1st_hvg_remove.xlsx')
hvg_keep = [gene for gene in hvg_list if gene not in remove_genes]
pd.DataFrame(hvg_keep).to_excel('./CD4T_1st_hvg_keep.xlsx')

is_keep = np.isin(adata.var_names, hvg_keep)
adata.var['highly_variable'] = is_keep
adata.var.highly_variable.value_counts()
adata_hvg = adata[:,adata.var["highly_variable"]==True].copy()

sc.pp.scale(adata_hvg,max_value=10)
sc.tl.pca(adata_hvg, n_comps = 50)
sc.pl.pca_variance_ratio(adata_hvg, log=True,n_pcs=50)
adata.uns['pca'] = adata_hvg.uns['pca'].copy()
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()

sc.pp.neighbors(adata, n_neighbors=30, n_pcs=25, use_rep='X_pca_harmony')
sc.tl.umap(adata,min_dist=0.3)
sc.tl.leiden(adata, resolution=1.5,n_iterations=2,key_added='leiden_r1.5_n2')

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['ExperimentalBatch'],
               legend_fontsize=8, legend_fontoutline=2)
    
with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['Date'],
                legend_fontsize=8, legend_fontoutline=2)
    
with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['leiden_r1.2_n2'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)
    
sc.pl.umap(adata, color=['CD3D','CD4','CD8A','CD8B',
                         'CD79A', 'HBB', 'PPBP','CD14', 'S100A8', 'MZB1',
                         'NCAM1', 'CD7', 'NCR1', 'IL7R',
                         'KLRC1', 'SLC4A10', 'FXYD2', 'HES1', 'NOTCH1', 'CD1A', 'SMPD3'])

# Isolation of residual CD8T
new_cluster_names = [
    'CD4T',
    'CD4T',
    'CD4T',
    'CD4T',
    'CD4T',
    'CD4T',#5
    'CD4T',
    'CD4T',
    'CD8T',
    'CD4T',
    'CD4T',#10
    
    'CD4T',#11
    'CD4T',
    'CD4T',
    'CD4T',
    'CD4T',#15
    'CD4T',
    'CD4T']
adata.obs['celltype_temp4']=adata.obs['leiden_r1.5_n2'].replace(to_replace = list(adata.obs['leiden_r1.5_n2'].cat.categories.values),
                                                         value = new_cluster_names)

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['celltype_temp4'],
               legend_fontsize=8, legend_fontoutline=2)

adata.obs.to_csv('./CD4T_temp1_metadata.csv')
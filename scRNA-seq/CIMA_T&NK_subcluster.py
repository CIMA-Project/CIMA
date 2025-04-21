import scanpy as sc
import anndata
import cupy as cp
import numpy as np
import pandas as pd
import time
import os
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

adata = sc.read_h5ad('cima/scRNA/CIMA_TNK_5.2M.h5ad')
adata.X = adata.layers['counts'].copy()

sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, batch_key='sample')
adata.var.highly_variable.value_counts()

mt_genes = list(adata.var_names[adata.var_names.str.match(r'^MT-')])
rp_genes = list(adata.var_names[adata.var_names.str.match(r'^RP[SL]')]) 
ncRNA_genes = list(adata.var_names[adata.var_names.str.match(r'^[A-Z][A-Z][0-9].*\.[0-9]')])
LINC_genes = list(adata.var_names[adata.var_names.str.match(r'(^LOC|LINC)[1-9]*')])
remove_genes = mt_genes + rp_genes + ncRNA_genes + LINC_genes

hvg_list = adata.var_names[adata.var.highly_variable].to_list()
pd.DataFrame(hvg_list).to_excel('CIMA/scRNA/TNK_hvg_list.xlsx')
hvg_remove = [gene for gene in hvg_list if gene in remove_genes]
pd.DataFrame(hvg_remove).to_excel('CIMA/scRNA/TNK_hvg_remove.xlsx')
hvg_keep = [gene for gene in hvg_list if gene not in remove_genes]
pd.DataFrame(hvg_keep).to_excel('CIMA/scRNA/TNK_hvg_keep.xlsx')

is_keep = np.isin(adata.var_names, hvg_keep)
adata.var['highly_variable'] = is_keep

adata_hvg = adata[:,adata.var["highly_variable"]==True].copy()
sc.pp.scale(adata_hvg,max_value=10)

sc.tl.pca(adata_hvg, n_comps = 50)
sc.pl.pca_variance_ratio(adata_hvg, log=True,n_pcs=50)

adata.uns['pca'] = adata_hvg.uns['pca'].copy()
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()

rsc.tl.harmony_integrate(adata,['Date'],max_iter_harmony = 30)
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=14,use_rep='X_pca_harmony')
sc.tl.umap(adata,min_dist=0.3)
sc.tl.leiden(adata, resolution=1.2,n_iterations=2,key_added='leiden_r1.2_n2')

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['ExperimentalBatch'],
               legend_fontsize=8, legend_fontoutline=2)
    
with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['Date'],
                legend_fontsize=8, legend_fontoutline=2)
    
with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['leiden_r1.2_n2'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)

sc.pl.umap(adata, color=['CD14', 'FCGR3A', 'FCN1',
                        'NCAM1', 'NCR1', 'KLRC1', 'KLRD1',
                        'IL7R','MKI67','CD3D', 'CD8A', 'CD4', 'CD79A', 'HBB', 'PPBP'])


# subclustering
adata_CD4T = adata[adata.obs['leiden_r1.2_n2'].isin(['0','3','1','14','18','17','13','4','6'])].copy()
adata_CD4T.write('CIMA/scRNA/CIMA_CD4T_doublet_2.4M.h5ad')

adata_CD8T = adata[~(adata.obs['leiden_r1.2_n2'].isin(['0','3','1','14','18','17','13','4','6']))].copy()
adata_CD8T.write('CIMA/scRNA/CIMA_CD48T_doublet_2.4M.h5ad')
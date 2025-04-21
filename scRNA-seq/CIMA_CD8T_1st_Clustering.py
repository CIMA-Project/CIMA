import os
os.environ["CUDA_VISIBLE_DEVICES"] = "1"
import scanpy as sc
import anndata
import cupy as cp
import numpy as np
import pandas as pd
import time
import rapids_singlecell as rsc
from rapids_singlecell.cunnData import cunnData

import warnings
warnings.filterwarnings("ignore")

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=150, facecolor='white') 

import rmm
rmm.reinitialize(
    managed_memory=True)
cp.cuda.set_allocator(rmm.rmm_cupy_allocator)

adata = sc.read_h5ad('CIMA/scRNA/CIMA_CD48T_doublet_2.4M.h5ad')
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
pd.DataFrame(hvg_list).to_excel('./CD8T_1st_hvg_list.xlsx')
hvg_remove = [gene for gene in hvg_list if gene in remove_genes]
pd.DataFrame(hvg_remove).to_excel('./CD8T_1st_hvg_remove.xlsx')
hvg_keep = [gene for gene in hvg_list if gene not in remove_genes]
pd.DataFrame(hvg_keep).to_excel('./CD8T_1st_hvg_keep.xlsx')

is_keep = np.isin(adata.var_names, hvg_keep)
adata.var['highly_variable'] = is_keep
adata.var.highly_variable.value_counts()
adata_hvg = adata[:,adata.var["highly_variable"]==True].copy()

sc.pp.scale(adata_hvg,max_value=10)
sc.tl.pca(adata_hvg, n_comps = 50)
sc.pl.pca_variance_ratio(adata_hvg, log=True,n_pcs=50)
adata.uns['pca'] = adata_hvg.uns['pca'].copy()
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()

rsc.tl.harmony_integrate(adata,['Date'],max_iter_harmony = 30)
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=24, use_rep='X_pca_harmony')
sc.tl.umap(adata,min_dist=0.3)
sc.tl.leiden(adata, resolution=1.5,n_iterations=2,key_added='leiden_r1.5_n2')

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['ExperimentalBatch'],
               legend_fontsize=8, legend_fontoutline=2)
    
with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['Date'],
               legend_fontsize=8, legend_fontoutline=2)
    
with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['leiden_r1.5_n2'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)
    
sc.pl.umap(adata, color=['CD3D','CD4','CD8A','CD8B',
                         'CD79A', 'HBB', 'PPBP','CD14', 'S100A8', 'MZB1',
                         'NCAM1', 'CD7', 'NCR1', 'IL7R',
                         'KLRC1', 'SLC4A10', 'FXYD2', 'HES1', 'NOTCH1', 'CD1A', 'SMPD3'])

# subset cluster 15 (Cycling T&NK) and reClustring to split Cycling T and Cycling NK
adata_cycling = adata[adata.obs['leiden_r1.5_n2']=='15'].copy()
adata_hvg = adata_cycling[:,adata_cycling.var["highly_variable"]==True].copy()

sc.pp.scale(adata_hvg,max_value=10)
sc.tl.pca(adata_hvg, n_comps = 50)
sc.pl.pca_variance_ratio(adata_hvg, log=True,n_pcs=50)
adata_cycling.uns['pca'] = adata_hvg.uns['pca'].copy()
adata_cycling.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()

rsc.tl.harmony_integrate(adata_cycling,['Date'],max_iter_harmony = 30)
sc.pp.neighbors(adata_cycling, n_neighbors=30, n_pcs=16, use_rep='X_pca_harmony')
sc.tl.umap(adata_cycling,min_dist=0.3)
sc.tl.leiden(adata_cycling, resolution=1.2,n_iterations=2,key_added='leiden_r1.2_n2')

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata_cycling, color=['ExperimentalBatch'],
               legend_fontsize=8, legend_fontoutline=2)
    
with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adadata_cyclingata, color=['Date'],
               legend_fontsize=8, legend_fontoutline=2)
    
with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata_cycling, color=['leiden_r1.2_n2'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)
    
sc.pl.dotplot(adata_cycling, ['CD3D','CD3G','CD4','CD8A','CD8B','CD7','NCAM1','NCR1','KLRC1','MKI67','TOP2A'], 'leiden_r1.2_n2', dendrogram=True)

sc.pl.umap(adata_cycling, color=['CD3D','CD4','CD8A','FCGR3A','MKI67','TOP2A','NCAM1','NCR1','CD7'])

new_cluster_names = [
    'Cycling_T',
    'Cycling_T',
    'Cycling_T',
    'Cycling_T',
    'Cycling_T',
    'Cycling_NK',#5
    'Cycling_T',
    'Cycling_T',
    'Cycling_T',
    'Cycling_NK',
    'Cycling_NK',#10
    'Cycling_T',
    'Cycling_T',
    'Cycling_T',
    'Cycling_T',
    'Cycling_T',#15
    'Cycling_NK',
    'Cycling_T',
    'Cycling_T']
adata_cycling.obs['celltype_temp3']=adata_cycling.obs['leiden_r1.2_n2'].replace(to_replace = list(adata_cycling.obs['leiden_r1.2_n2'].cat.categories.values),
                                                         value = new_cluster_names)

adata_cycling.obs.to_csv('CD8T_temp3_metadata.csv')


# Isolation of residual CD4T
adata_0_5_13 = adata[adata.obs['leiden_r1.5_n2'].isin(['5','0','13'])].copy()
adata_hvg = adata_0_5_13[:,adata_0_5_13.var["highly_variable"]==True].copy()

sc.pp.scale(adata_hvg,max_value=10)
sc.tl.pca(adata_hvg, n_comps = 50)
sc.pl.pca_variance_ratio(adata_hvg, log=True,n_pcs=50)
adata_0_5_13.uns['pca'] = adata_hvg.uns['pca'].copy()
adata_0_5_13.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()

rsc.tl.harmony_integrate(adata_0_5_13,['Date'],max_iter_harmony = 30)
sc.pp.neighbors(adata_0_5_13, n_neighbors=30, n_pcs=25, use_rep='X_pca_harmony')
sc.tl.umap(adata_0_5_13,min_dist=0.3)
sc.tl.leiden(adata_0_5_13, resolution=1.2,n_iterations=2,key_added='leiden_r1.2_n2')

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata_0_5_13, color=['ExperimentalBatch'],
               legend_fontsize=8, legend_fontoutline=2)
    
with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata_0_5_13, color=['Date'],
               legend_fontsize=8, legend_fontoutline=2)
    
with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata_0_5_13, color=['leiden_r1.2_n2'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)
    
sc.pl.umap(adata_0_5_13, color=['CD3D','CD4','CD8A','FCGR3A','MKI67'])
sc.pl.umap(adata_0_5_13[adata_0_5_13.obs['leiden_r1.2_n2'].isin(['4','0','5','3','8','6'])], color='CD4')

adata_temp = adata_0_5_13[adata_0_5_13.obs['leiden_r1.2_n2'].isin(['4','0','5','3','8','6'])].copy()
adata_hvg = adata_temp[:,adata_temp.var["highly_variable"]==True].copy()

sc.pp.scale(adata_hvg,max_value=10)
sc.tl.pca(adata_hvg, n_comps = 50)
sc.pl.pca_variance_ratio(adata_hvg, log=True,n_pcs=50)
adata_temp.uns['pca'] = adata_hvg.uns['pca'].copy()
adata_temp.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()

rsc.tl.harmony_integrate(adata_temp,['Date'],max_iter_harmony = 30)
sc.pp.neighbors(adata_temp, n_neighbors=30, n_pcs=13, use_rep='X_pca_harmony')
sc.tl.umap(adata_temp,min_dist=0.3)
sc.tl.leiden(adata_temp, resolution=1.5,n_iterations=2,key_added='leiden_r1.5_n2')

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata_temp, color=['ExperimentalBatch'],
               legend_fontsize=8, legend_fontoutline=2)
    
with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata_temp, color=['Date'],
               legend_fontsize=8, legend_fontoutline=2)
    
with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata_temp, color=['leiden_r1.2_n2'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)
    
sc.pl.umap(adata_temp, color=['CD3D','CD4','CD8A','FCGR3A','MKI67'])

adata_temp = adata_temp[adata_temp.obs['leiden_r1.5_n2'].isin(['9','2','4','0'])].copy()
adata_hvg = adata_temp[:,adata_temp.var["highly_variable"]==True].copy()

sc.pp.scale(adata_hvg,max_value=10)
sc.tl.pca(adata_hvg, n_comps = 50)
sc.pl.pca_variance_ratio(adata_hvg, log=True,n_pcs=50)
adata_temp.uns['pca'] = adata_hvg.uns['pca'].copy()
adata_temp.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()

rsc.tl.harmony_integrate(adata_temp,['Date'],max_iter_harmony = 30)
sc.pp.neighbors(adata_temp, n_neighbors=30, n_pcs=15, use_rep='X_pca_harmony')
sc.tl.umap(adata_temp,min_dist=0.3)
sc.tl.leiden(adata_temp, resolution=1.2,n_iterations=2,key_added='leiden_r1.2_n2')

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata_temp, color=['leiden_r1.2_n2'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)

sc.pl.umap(adata_temp, color=['CD3D','CD4','CD8A','FCGR3A','MKI67'])

new_cluster_names = [
    'CD8T',
    'CD4T',
    'CD4T',
    'CD4T',
    'CD4T',
    'CD4T',#5
    'CD4T',
    'CD4T']
adata_temp.obs['celltype_temp1']=adata_temp.obs['leiden_r1.2_n2'].replace(to_replace = list(adata_temp.obs['leiden_r1.2_n2'].cat.categories.values),
                                                         value = new_cluster_names)

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata_temp, color=['celltype_temp1'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)
    
adata_temp.obs.to_csv('CD8T_temp1_metadata.csv')


# Isolation of residual NK
adata_9 = adata[adata.obs['leiden_r1.5_n2'].isin(['9'])].copy()
adata_hvg = adata_9[:,adata_9.var["highly_variable"]==True].copy()

sc.pp.scale(adata_hvg,max_value=10)
sc.tl.pca(adata_hvg, n_comps = 50)
sc.pl.pca_variance_ratio(adata_hvg, log=True,n_pcs=50)
adata_9.uns['pca'] = adata_hvg.uns['pca'].copy()
adata_9.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()

rsc.tl.harmony_integrate(adata_9,['Date'],max_iter_harmony = 30)
sc.pp.neighbors(adata_9, n_neighbors=30, n_pcs=18, use_rep='X_pca_harmony')
sc.tl.umap(adata_9,min_dist=0.3)
sc.tl.leiden(adata_9, resolution=1.2,n_iterations=2,key_added='leiden_r1.2_n2')

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata_9, color=['leiden_r1.2_n2'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)
    
sc.pl.umap(adata_9,color=['CD3D', 'CD3G', 'CD4', 'CD8A', 'CD8B', 'CD7', 'NCAM1', 'FCGR3A',  'KLRC1', 'KLRD1', 'KLRF1', 'NKG7', 'GNLY'])

new_cluster_names = [
    'CD8T',
    'CD8T',
    'CD8T',
    'NK',
    'NK',
    'NK',#5
    'NK',
    'CD8T',
    'NK',
    'CD8T',
    'NK',#10
    'NK',
    'CD8T']
adata_9.obs['celltype_temp2']=adata_9.obs['leiden_r1.2_n2'].replace(to_replace = list(adata_9.obs['leiden_r1.2_n2'].cat.categories.values),
                                                         value = new_cluster_names)

adata_9.obs.to_csv('CD8T_temp2_metadata.csv')
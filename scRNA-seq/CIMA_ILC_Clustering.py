import os
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
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

adata = sc.read_h5ad('/media/AnalysisFastDisk/NatualCohort_NK_0.5M.h5ad')
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
pd.DataFrame(hvg_list).to_excel('./NK_1st_hvg_list.xlsx')
hvg_remove = [gene for gene in hvg_list if gene in remove_genes]
pd.DataFrame(hvg_remove).to_excel('./NK_1st_hvg_remove.xlsx')
hvg_keep = [gene for gene in hvg_list if gene not in remove_genes]
pd.DataFrame(hvg_keep).to_excel('./NK_1st_hvg_keep.xlsx')

is_keep = np.isin(adata.var_names, hvg_keep)
adata.var['highly_variable'] = is_keep
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
    
NK_marker_genes = ['CD3D', 'CD3G', 'CD4', 'CD8A', 'CD8B', 'CD7', 'NCAM1', 'FCGR3A',  'KLRC1', 'KLRD1', 'KLRF1',
                 'NKG7', 'GNLY', 'CD160', 'XCL1', 'XCL2', 'GZMA', 'GZMK', 'AXAN1', 'AXAN2', 'ACTB', 'S100A4', 'NR4A2',
                 'DUSP1', 'FOSB', 'SIGLEC7', 'IFNAR2', 'CCL3', 'CCL4', 'ISG15', 'KLRB1', 'LILRB1', 'KIR2DL1', 'B3GAT1', 'CX3CR1',
                 'HAVCR2', 'ZEB2', 'WDR74', 'IL7R', 'KIT', 'TBX21', 'GATA3', 'KLRG1', 'IL7RB', 'PCDH9', 'LST1', 'NCR1', 'NCR2']
NK_marker_genes = list(filter(lambda x: x in adata.var_names, NK_marker_genes))
len(NK_marker_genes)
sc.pl.umap(adata, color=NK_marker_genes)

# remove Cells contaminated with red blood cells
adata = adata[~(adata.obs['leiden_r1.5_n2'].isin(['9']))].copy()
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
pd.DataFrame(hvg_list).to_excel('./NK_2nd_hvg_list.xlsx')
hvg_remove = [gene for gene in hvg_list if gene in remove_genes]
pd.DataFrame(hvg_remove).to_excel('./NK_2nd_hvg_remove.xlsx')
hvg_keep = [gene for gene in hvg_list if gene not in remove_genes]
pd.DataFrame(hvg_keep).to_excel('./NK_2nd_hvg_keep.xlsx')

is_keep = np.isin(adata.var_names, hvg_keep)
adata.var['highly_variable'] = is_keep
adata_hvg = adata[:,adata.var["highly_variable"]==True].copy()

sc.pp.scale(adata_hvg,max_value=10)
sc.tl.pca(adata_hvg, n_comps = 50)
sc.pl.pca_variance_ratio(adata_hvg, log=True,n_pcs=50)
adata.uns['pca'] = adata_hvg.uns['pca'].copy()
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()

rsc.tl.harmony_integrate(adata,['Date'],max_iter_harmony = 30)
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=22, use_rep='X_pca_harmony')
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
    
clusters = list(adata.obs['leiden_r1.5_n2'].cat.categories.values)
for cluster in clusters:
    sc.pl.umap(adata,color=['leiden_r1.5_n2'],groups=[cluster],save='_NK_2nd_cluster_'+cluster+'.png')

NK_marker_genes = ['CD3D', 'CD3G', 'CD4', 'CD8A', 'CD8B', 'CD7', 'NCAM1', 'FCGR3A',  'KLRC1', 'KLRD1', 'KLRF1',
                 'NKG7', 'GNLY', 'CD160', 'XCL1', 'XCL2', 'GZMA', 'GZMK', 'AXAN1', 'AXAN2', 'ACTB', 'S100A4', 'NR4A2',
                 'DUSP1', 'FOSB', 'SIGLEC7', 'IFNAR2', 'CCL3', 'CCL4', 'ISG15', 'KLRB1', 'LILRB1', 'KIR2DL1', 'B3GAT1', 'CX3CR1',
                 'HAVCR2', 'ZEB2', 'WDR74', 'IL7R', 'KIT', 'TBX21', 'GATA3', 'KLRG1', 'IL7RB', 'PCDH9', 'LST1', 'NCR1', 'NCR2']
NK_marker_genes = list(filter(lambda x: x in adata.var_names, NK_marker_genes))
len(NK_marker_genes)
sc.pl.umap(adata, color=NK_marker_genes)

adata.write('CIMA/scRNA/CIMA_Cluster2nd.h5ad')


# Reclustering to obtain accurated clusters, based on marker gene exp.
adata_temp1 = adata[adata.obs['leiden_r1.5_n2'].isin(['2','6','7','10'])].copy()
adata_hvg = adata_temp1[:,adata_temp1.var["highly_variable"]==True].copy()

sc.pp.scale(adata_hvg,max_value=10)
sc.tl.pca(adata_hvg, n_comps = 50)
sc.pl.pca_variance_ratio(adata_hvg, log=True,n_pcs=50)
adata_temp1.uns['pca'] = adata_hvg.uns['pca'].copy()
adata_temp1.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()

rsc.tl.harmony_integrate(adata_temp1,['Date'],max_iter_harmony = 30)
sc.pp.neighbors(adata_temp1, n_neighbors=30, n_pcs=13, use_rep='X_pca_harmony')
sc.tl.umap(adata_temp1,min_dist=0.3)
sc.tl.leiden(adata_temp1, resolution=1.5,n_iterations=2,key_added='leiden_r1.5_n2')

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata_temp1, color=['ExperimentalBatch'],
               legend_fontsize=8, legend_fontoutline=2)
    
with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata_temp1, color=['Date'],
                legend_fontsize=8, legend_fontoutline=2)
    
with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata_temp1, color=['leiden_r1.5_n2'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)
    
sc.pl.umap(adata_temp1, color=NK_marker_genes)

new_cluster_names = [
    'NK_ZEB2',
    'NK_ZEB2',
    'NK_ZEB2',
    'NK_GZMK',
    'ILC_NCAM1',
    'NK_GZMK',#5
    'NK_GZMK',
    'NK_GZMK',
    'ILC_NCAM1',
    'NK_GZMK',
    'ILC_NCAM1',#10
    
    'ILC_KIT',#11
    'NK_ZEB2']
adata_temp1.obs['celltype_NK_temp']=adata_temp1.obs['leiden_r1.5_n2'].replace(to_replace = list(adata_temp1.obs['leiden_r1.5_n2'].cat.categories.values),
                                                         value = new_cluster_names)

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata_temp1, color=['celltype_NK_temp'],legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)
    
adata_temp1.write('CIMA/scRNA/CIMA_NK_temp1.h5ad')

adata.obs['celltype_NK_1st'] = adata.obs['leiden_r1.5_n2']
adata.obs['celltype_NK_1st'] = adata.obs.index.map(adata_temp1.obs['celltype_NK_temp'].to_dict())
adata.obs['celltype_NK_1st'].fillna(adata.obs['leiden_r1.5_n2'],inplace=True)

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['celltype_NK_1st'],legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)
    
adata.write('CIMA/scRNA/CIMA_ILC_Annotation.h5ad')
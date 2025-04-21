import scanpy as sc
import os
import anndata as ad
import numpy as np
import pandas as pd
import seaborn as sb
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import rapids_singlecell as rsc
import cupy as cp

import warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=150, facecolor='white') 

import rmm
rmm.reinitialize(
    managed_memory=True)
cp.cuda.set_allocator(rmm.rmm_cupy_allocator)

adata = sc.read_h5ad('CIMA/scRNA/CIMA_B_cells_0.6M.h5ad')
adata.obs = adata.obs.drop(['leiden_r1.5_n2', 'leiden_r2_n2'],axis = 1)
adata.var = adata.var.drop(['highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm'],axis = 1)
del adata.uns
del adata.obsm
del adata.obsp
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
pd.DataFrame(hvg_list).to_excel('CIMA/scRNA/B_hvg_list.xlsx')
hvg_remove = [gene for gene in hvg_list if gene in remove_genes]
pd.DataFrame(hvg_remove).to_excel('CIMA/scRNA/B_hvg_remove.xlsx')
hvg_keep = [gene for gene in hvg_list if gene not in remove_genes]
pd.DataFrame(hvg_keep).to_excel('CIMA/scRNA/B_hvg_keep.xlsx')

is_keep = np.isin(adata.var_names, hvg_keep)
adata.var['highly_variable'] = is_keep

adata_hvg = adata[:,adata.var["highly_variable"]==True].copy()
sc.pp.scale(adata_hvg,max_value=10)

rsc.tl.pca(adata_hvg, n_comps = 50)
sc.pl.pca_variance_ratio(adata_hvg, log=True,n_pcs=50)
sc.pl.pca_loadings(adata_hvg, components=range(1,30))

adata.uns['pca'] = adata_hvg.uns['pca'].copy()
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()

rsc.tl.harmony_integrate(adata,['Date'],max_iter_harmony = 30)
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=15, use_rep='X_pca_harmony')
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
    sc.pl.umap(adata,color=['leiden_r1.5_n2'],groups=[cluster],save='_B_1st_cluster_'+cluster+'.pdf')

sc.pl.umap(adata, color=['MKI67','CD3D','CD19','XBP1','CD79A',
                         'IGHA1','IGHG1'])

new_cluster_names = [
    'B cells',
    'B cells',
    'B cells',
    'B cells',
    'B cells',
    'B cells',#5
    'B cells',
    'B cells',
    'B cells',
    'B cells',
    'B cells',#10
    
    'B cells',#11
    'B cells',
    'B cells',
    'doublet',
    'Plasma',#15
    'Cycling T',
    'doublet',
    'Cycling T',
    'Plasma',
    'Plasma']
adata.obs['celltype_B_1st']=adata.obs['leiden_r1.5_n2'].replace(to_replace = list(adata.obs['leiden_r1.5_n2'].cat.categories.values),
                                                         value = new_cluster_names)

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['celltype_B_1st'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['celltype_1st'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)
    
adata.obs.to_csv('CIMA/scRNA/B_Annotation1st_Metadata.csv')
adata.write('CIMA/scRNA/CIMA_B_Annotation1st.h5ad')

# remove doublet and Cycling T
adata = adata[adata.obs['celltype_B_1st'].isin(['B cells','Plasma'])].copy()
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
pd.DataFrame(hvg_list).to_excel('CIMA/scRNA/B_singlet_hvg_list.xlsx')
hvg_remove = [gene for gene in hvg_list if gene in remove_genes]
pd.DataFrame(hvg_remove).to_excel('CIMA/scRNA/B_singlet_hvg_remove.xlsx')
hvg_keep = [gene for gene in hvg_list if gene not in remove_genes]
pd.DataFrame(hvg_keep).to_excel('CIMA/scRNA/B_singlet_hvg_keep.xlsx')

is_keep = np.isin(adata.var_names, hvg_keep)
adata.var['highly_variable'] = is_keep

adata_hvg = adata[:,adata.var["highly_variable"]==True].copy()
sc.pp.scale(adata_hvg,max_value=10)

rsc.tl.pca(adata_hvg, n_comps = 50)
sc.pl.pca_variance_ratio(adata_hvg, log=True,n_pcs=50)
sc.pl.pca_loadings(adata_hvg, components=range(1,30))

adata.uns['pca'] = adata_hvg.uns['pca'].copy()
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()

rsc.tl.harmony_integrate(adata,['Date'],max_iter_harmony = 30)
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=13, use_rep='X_pca_harmony')
sc.tl.umap(adata,min_dist=0.5)
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
    
sc.pl.umap(adata, color=['MKI67','CD3D','CD19','XBP1','CD79A',
                         'IGHA1','IGHG1'])

clusters = list(adata.obs['leiden_r1.5_n2'].cat.categories.values)
for cluster in clusters:
    sc.pl.umap(adata,color=['leiden_r1.5_n2'],groups=[cluster],save='_B_2nd_cluster_'+cluster+'.pdf')

B_marker_genes = ['MS4A1', 'SDC1', 'NEIL1', 'TCL1A', 'AICDA', 'BCL6', 'CD24', 'MYO1C',
                  'IGHM', 'IGHD', 'IL4R', 'CCR7', 'CD27', 'CXCR5', 'CR2', 'ITGAX','TBX21',
                  'FCRL2', 'FCRL3', 'FCRL4', 'FCRL5', 'CD1C', 'TNFSF13B', 'CD86', 'ZEB2', 
                  'FGR', 'CIB1', 'CD22', 'CD69', 'FCER2', 'CD38', 'MCM4', 'PCNA', 'MKI67',
                  'CERCAM', 'ITGA8', 'IGHG1', 'IGHA1', 'HLA-DRB5']
B_marker_genes = list(filter(lambda x: x in adata.var_names, B_marker_genes))
len(B_marker_genes)
sc.pl.umap(adata, color=B_marker_genes)

genes_keep = [gene for gene in adata.var_names if gene not in remove_genes]
len(genes_keep)
adata_for_deg = adata[:,genes_keep].copy()

rsc.tl.rank_genes_groups_logreg(adata_for_deg, groupby="leiden_r1.5_n2", use_raw=False)
sc.pl.rank_genes_groups(adata_for_deg, n_genes=20)

# COSG
import cosg as cosg
import importlib
importlib.reload(cosg)

cosg.cosg(adata_for_deg,
    key_added='cosg',
    mu=1,
    n_genes_user=50,
    groupby='leiden_r1.5_n2')

colnames = ['names', 'scores']
test = [pd.DataFrame(adata_for_deg.uns["cosg"][c]) for c in colnames]
test = pd.concat(test, axis=1, names=[None, 'group'], keys=colnames)

markers = {}
cats = adata_for_deg.obs['leiden_r1.5_n2'].cat.categories
for i, c in enumerate(cats):
    cell_type_df = test.loc[:, 'names'][str(c)]
    scores_df = test.loc[:, 'scores'][str(c)]
    markers[str(c)] = cell_type_df.values.tolist()[:10]

sc.pl.dotplot(adata_for_deg, 
              var_names = markers, 
              groupby='leiden_r1.5_n2',
              cmap='Spectral_r',
              standard_scale='var',save='_B_singlet_leiden_r1.5_n2.pdf') 

marker_gene=pd.DataFrame(adata_for_deg.uns['cosg']['names'])
marker_gene.to_excel('CIMA/scRNA/B_singlet_leiden_r1.5_n2_marker_gene.xlsx')

marker_gene_scores=pd.DataFrame(adata_for_deg.uns['cosg']['scores'])
marker_gene_scores.to_excel('CIMA/scRNA/B_singlet_leiden_r1.5_n2_marker_gene_scores.xlsx')

adata.write('CIMA/scRNA/CIMA_B_Cluster2nd_0.5M.h5ad')

# Annotation
new_cluster_names = [
    'Bn_IGHM',
    'Switched_resting_Bm_CD27',
    'Bn_IL4R',
    'Bn_IL4R',
    'Transitional_B_CD9',
    'Activated_Bn_FOS',#5
    'Unswitched_Bm_CD1C',
    'Switched_Bm_IGHE',
    'Switched_activated_Bm_CD86',
    'Pre-switch_Bm_JAM3',
    'Atypical_Bm_ITGAX',#10
    
    'Switched_Bm_IGHA1',#11
    'Bn_IL6',
    'Bn_IFIT3',
    'Plasma_IGHA1',
    'Transitional_B_MME',#15
    'Bn_IGHD',
    'Unswitched_Bm_IL6',
    'Plasma_IGHG1',
    'Plasmablast_MKI67',
    'Pre-switch_Bm_IFIT3']
adata.obs['celltype_B_2nd']=adata.obs['leiden_r1.5_n2'].replace(to_replace = list(adata.obs['leiden_r1.5_n2'].cat.categories.values),
                                                         value = new_cluster_names)

with rc_context({'figure.figsize': (15, 15)}):
    sc.pl.umap(adata, color=['celltype_B_2nd'],palette=B_colors,
               legend_fontsize=8, legend_fontoutline=2,size=5)
    
adata.write('CIMA/scRNA/CIMA_B_Annotation_2nd.h5ad')
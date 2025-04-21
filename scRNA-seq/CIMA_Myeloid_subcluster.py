import scanpy as sc
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "1"
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

adata = sc.read_h5ad('CIMA/scRNA/CIMA_myeloid_0.9M.h5ad')
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
pd.DataFrame(hvg_list).to_excel('CIMA/scRNA/Myeloid_hvg_list.xlsx')
hvg_remove = [gene for gene in hvg_list if gene in remove_genes]
pd.DataFrame(hvg_remove).to_excel('CIMA/scRNA/Myeloid_hvg_remove.xlsx')
hvg_keep = [gene for gene in hvg_list if gene not in remove_genes]
pd.DataFrame(hvg_keep).to_excel('CIMA/scRNA/Myeloid_hvg_keep.xlsx')

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
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=21, use_rep='X_pca_harmony')
sc.tl.umap(adata,min_dist=0.3)
sc.tl.leiden(adata, resolution=1.8,n_iterations=2,key_added='leiden_r1.8_n2')

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['ExperimentalBatch'],
               legend_fontsize=8, legend_fontoutline=2)
    
with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['Date'],
               legend_fontsize=8, legend_fontoutline=2)
    
with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['leiden_r1.8_n2'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)
    
clusters = list(adata.obs['leiden_r1.8_n2'].cat.categories.values)
for cluster in clusters:
    sc.pl.umap(adata,color=['leiden_r1.8_n2'],groups=[cluster],save='_myeloid_1st_cluster_'+cluster+'.pdf')

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['celltype_1st'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)
    
sc.pl.umap(adata, color=['CD3D', 'CD8A', 'CD4', 'CD79A', 'HBB', 'PPBP'])

myeloid_marker_genes = ['CD74', 'HLA-DPB1', 'CLEC9A', 'BATF3', 'IDO1', 'CD1C', 'LAMP3',
                        'CLEC10A', 'PKIB', 'ITGAX', 'FCGR3A', 'AXL', 'SIGLEC6', 'FCN1', 
                        'CD14', 'CD86', 'CRISPLD2', 'CLU', 'MXD1', 'CXCR1', 'CD83', 'CPA3',
                        'TPSAB1', 'KIT', 'CYTL1', 'SOX4', 'ITGAM', 'CEACAM8', 'CD81', 'SPN', 'FUT4', 'MME']
len(myeloid_marker_genes)
myeloid_marker_genes = list(filter(lambda x: x in adata.var_names, myeloid_marker_genes))
len(myeloid_marker_genes)
sc.pl.umap(adata, color=myeloid_marker_genes)

genes_keep = [gene for gene in adata.var_names if gene not in remove_genes]
len(genes_keep)
adata_for_deg = adata[:,genes_keep].copy()
rsc.tl.rank_genes_groups_logreg(adata_for_deg, groupby="leiden_r1.8_n2", use_raw=False)
sc.pl.rank_genes_groups(adata_for_deg, n_genes=20)

# COSG
import cosg as cosg
import importlib
importlib.reload(cosg)

cosg.cosg(adata_for_deg,
    key_added='cosg',
    mu=1,
    n_genes_user=50,
    groupby='leiden_r1.8_n2')

colnames = ['names', 'scores']
test = [pd.DataFrame(adata_for_deg.uns["cosg"][c]) for c in colnames]
test = pd.concat(test, axis=1, names=[None, 'group'], keys=colnames)

markers = {}
cats = adata_for_deg.obs['leiden_r1.8_n2'].cat.categories
for i, c in enumerate(cats):
    cell_type_df = test.loc[:, 'names'][str(c)]
    scores_df = test.loc[:, 'scores'][str(c)]
    markers[str(c)] = cell_type_df.values.tolist()[:10]

sc.pl.dotplot(adata_for_deg, 
              var_names = markers, 
              groupby='leiden_r1.8_n2',
              cmap='Spectral_r',
              standard_scale='var',save='_Myeloid1st_dotplot_leiden_r1.8_n2.pdf') 

marker_gene=pd.DataFrame(adata_for_deg.uns['cosg']['names'])
marker_gene.to_excel('CIMA/scRNA/Myeloid1st_dotplot_leiden_r1.8_n2_marker_gene.xlsx')

marker_gene_scores=pd.DataFrame(adata_for_deg.uns['cosg']['scores'])
marker_gene_scores.to_excel('CIMA/scRNA/Myeloid1st_dotplot_leiden_r1.8_n2_marker_gene_scores.xlsx')

adata.write('CIMA/scRNA/NatualCohort_myeloid_Cluster1st_0.9M.h5ad')


# remove doublet round 1 and 2nd clustering
adata = adata[~(adata.obs['leiden_r1.8_n2'].isin(['1','11','19']))].copy()
adata.X = adata.layers['counts'].copy()

sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, batch_key="sample")
adata.var.highly_variable.value_counts()

mt_genes = list(adata.var_names[adata.var_names.str.match(r'^MT-')])
rp_genes = list(adata.var_names[adata.var_names.str.match(r'^RP[SL]')]) 
ncRNA_genes = list(adata.var_names[adata.var_names.str.match(r'^[A-Z][A-Z][0-9].*\.[0-9]')])
LINC_genes = list(adata.var_names[adata.var_names.str.match(r'(^LOC|LINC)[1-9]*')])
remove_genes = mt_genes + rp_genes + ncRNA_genes + LINC_genes

hvg_list = adata.var_names[adata.var.highly_variable].to_list()
pd.DataFrame(hvg_list).to_excel('CIMA/scRNA/Myeloid_2nd_hvg_list.xlsx')
hvg_remove = [gene for gene in hvg_list if gene in remove_genes]
pd.DataFrame(hvg_remove).to_excel('CIMA/scRNA/Myeloid_2nd_hvg_remove.xlsx')
hvg_keep = [gene for gene in hvg_list if gene not in remove_genes]
pd.DataFrame(hvg_keep).to_excel('CIMA/scRNA/Myeloid_2nd_hvg_keep.xlsx')

is_keep = np.isin(adata.var_names, hvg_keep)
adata.var['highly_variable'] = is_keep

adata_hvg = adata[:,adata.var["highly_variable"]==True].copy()
sc.pp.scale(adata_hvg,max_value=10)

rsc.tl.pca(adata_hvg, n_comps = 50)
sc.pl.pca_variance_ratio(adata_hvg, log=True,n_pcs=50)

adata.uns['pca'] = adata_hvg.uns['pca'].copy()
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()

rsc.tl.harmony_integrate(adata,'Date',max_iter_harmony = 30)
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=16, use_rep='X_pca_harmony')
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
    
clusters = list(adata.obs['leiden_r1.5_n2'].cat.categories.values)
for cluster in clusters:
    sc.pl.umap(adata,color=['leiden_r1.5_n2'],groups=[cluster],save='_myeloid_2nd_cluster_'+cluster+'.png')

sc.pl.umap(adata, color=['CD3D', 'CD8A', 'CD4', 'CD79A', 'HBB', 'PPBP'])

myeloid_marker_genes = ['CD74', 'HLA-DPB1', 'CLEC9A', 'BATF3', 'IDO1', 'CD1C', 'LAMP3',
                        'CLEC10A', 'PKIB', 'ITGAX', 'FCGR3A', 'AXL', 'SIGLEC6', 'FCN1', 
                        'CD14', 'CD86', 'CRISPLD2', 'CLU', 'MXD1', 'CXCR1', 'CD83', 'CPA3',
                        'TPSAB1', 'KIT', 'CYTL1', 'SOX4', 'ITGAM', 'CEACAM8', 'CD81', 'SPN', 'FUT4', 'MME']
len(myeloid_marker_genes)
myeloid_marker_genes = list(filter(lambda x: x in adata.var_names, myeloid_marker_genes))
len(myeloid_marker_genes)
sc.pl.umap(adata, color=myeloid_marker_genes)

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
              standard_scale='var',save='_Myeloid2nd_dotplot_leiden_r1.5_n2.pdf') 

marker_gene=pd.DataFrame(adata_for_deg.uns['cosg']['names'])
marker_gene.to_excel('CIMA/scRNA/Myeloid2nd_dotplot_leiden_r1.8_n2_marker_gene.xlsx')

marker_gene_scores=pd.DataFrame(adata_for_deg.uns['cosg']['scores'])
marker_gene_scores.to_excel('CIMA/scRNA/Myeloid2nd_dotplot_leiden_r1.8_n2_marker_gene_scores.xlsx')

adata.write('CIMA/scRNA/CIMA_myeloid_Cluster2nd_0.7M.h5ad')

# remove doublet round 2 and 3rd clustering
# 19: T&NK cell  22: doublet
adata = adata[~(adata.obs['leiden_r1.5_n2'].isin(['19','22']))].copy()
adata.X = adata.layers['counts'].copy()

sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, batch_key="sample")
adata.var.highly_variable.value_counts()

mt_genes = list(adata.var_names[adata.var_names.str.match(r'^MT-')])
rp_genes = list(adata.var_names[adata.var_names.str.match(r'^RP[SL]')]) 
ncRNA_genes = list(adata.var_names[adata.var_names.str.match(r'^[A-Z][A-Z][0-9].*\.[0-9]')])
LINC_genes = list(adata.var_names[adata.var_names.str.match(r'(^LOC|LINC)[1-9]*')])
remove_genes = mt_genes + rp_genes + ncRNA_genes + LINC_genes

hvg_list = adata.var_names[adata.var.highly_variable].to_list()
pd.DataFrame(hvg_list).to_excel('CIMA/scRNA/Myeloid_3rd_hvg_list.xlsx')
hvg_remove = [gene for gene in hvg_list if gene in remove_genes]
pd.DataFrame(hvg_remove).to_excel('CIMA/scRNA/Myeloid_3rd_hvg_remove.xlsx')
hvg_keep = [gene for gene in hvg_list if gene not in remove_genes]
pd.DataFrame(hvg_keep).to_excel('CIMA/scRNA/Myeloid_3rd_hvg_keep.xlsx')

is_keep = np.isin(adata.var_names, hvg_keep)
adata.var['highly_variable'] = is_keep

adata_hvg = adata[:,adata.var["highly_variable"]==True].copy()
sc.pp.scale(adata_hvg,max_value=10)

rsc.tl.pca(adata_hvg, n_comps = 50)
sc.pl.pca_variance_ratio(adata_hvg, log=True,n_pcs=50)

adata.uns['pca'] = adata_hvg.uns['pca'].copy()
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()

rsc.tl.harmony_integrate(adata,'Date',max_iter_harmony = 30)
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=20, use_rep='X_pca_harmony')
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
    
clusters = list(adata.obs['leiden_r1.5_n2'].cat.categories.values)
for cluster in clusters:
    sc.pl.umap(adata,color=['leiden_r1.5_n2'],groups=[cluster],save='_myeloid_3rd_cluster_'+cluster+'.png')

sc.pl.umap(adata, color=['CD3D', 'CD8A', 'CD4', 'CD79A', 'HBB', 'PPBP'])

myeloid_marker_genes = ['CD74', 'HLA-DPB1', 'CLEC9A', 'BATF3', 'IDO1', 'CD1C', 'LAMP3',
                        'CLEC10A', 'PKIB', 'ITGAX', 'FCGR3A', 'AXL', 'SIGLEC6', 'FCN1', 
                        'CD14', 'CD86', 'CRISPLD2', 'CLU', 'MXD1', 'CXCR1', 'CD83', 'CPA3',
                        'TPSAB1', 'KIT', 'CYTL1', 'SOX4', 'ITGAM', 'CEACAM8', 'CD81', 'SPN', 'FUT4', 'MME']
len(myeloid_marker_genes)
myeloid_marker_genes = list(filter(lambda x: x in adata.var_names, myeloid_marker_genes))
len(myeloid_marker_genes)
sc.pl.umap(adata, color=myeloid_marker_genes)

adata.write('CIMA/scRNA/CIMA_myeloid_Cluster3rd_0.7M.h5ad')

# subclustering 13 and 20 clusters in the 3rd cluster for more accurate annotation
adata_sub = adata[adata.obs['leiden_r1.5_n2'].isin(['13','20'])].copy()
adata_hvg = adata_sub[:,adata_sub.var["highly_variable"]==True].copy()

sc.pp.scale(adata_hvg,max_value=10)
sc.tl.pca(adata_hvg, n_comps = 50)
sc.pl.pca_variance_ratio(adata_hvg, log=True,n_pcs=50)

adata_sub.uns['pca'] = adata_hvg.uns['pca'].copy()
adata_sub.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()

rsc.tl.harmony_integrate(adata_sub,'Date',max_iter_harmony = 30)
sc.pp.neighbors(adata_sub, n_neighbors=30, n_pcs=13, use_rep='X_pca_harmony')
sc.tl.umap(adata_sub,min_dist=0.5)
sc.tl.leiden(adata_sub, resolution=0.5,key_added='leiden_r0.5')

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata_sub, color=['leiden_r0.5'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)
    
sc.pl.umap(adata_sub, color=['CD3D', 'CD8A', 'CD4', 'CD79A', 'HBB', 'PPBP'])
sc.pl.umap(adata_sub, color=['LILRA4','GZMB','IRF4', 'NCAM1', 'HLA-DPA1', 'LAD1', 'FSCN1', 'CCL22'])
myeloid_marker_genes = ['CD74', 'HLA-DPB1', 'CLEC9A', 'BATF3', 'IDO1', 'CD1C', 'LAMP3',
                        'CLEC10A', 'PKIB', 'ITGAX', 'FCGR3A', 'AXL', 'SIGLEC6', 'FCN1', 
                        'CD14', 'CD86', 'CRISPLD2', 'CLU', 'MXD1', 'CXCR1', 'CD83', 'CPA3',
                        'TPSAB1', 'KIT', 'CYTL1', 'SOX4', 'ITGAM', 'CEACAM8', 'CD81', 'SPN', 'FUT4', 'MME']
len(myeloid_marker_genes)
myeloid_marker_genes = list(filter(lambda x: x in adata_sub.var_names, myeloid_marker_genes))
len(myeloid_marker_genes)
sc.pl.umap(adata_sub, color=myeloid_marker_genes)

new_cluster_names = [
    'pDC_LILRA4',
    'pDC_LILRA4',
    'pDC_LILRA4',
    'pDC',
    'doublet',
    'pDC_LILRA4',#5
    'AS_DC',
    'doublet',
    'pDC_FCN1']
adata_sub.obs['leiden_r0.5']=adata_sub.obs['leiden_r0.5'].replace(to_replace = list(adata_sub.obs['leiden_r0.5'].cat.categories.values),
                                                         value = new_cluster_names)

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata_sub, color=['leiden_r0.5'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)
    
adata.obs['celltype_sub'] = adata.obs['leiden_r1.5_n2']
adata.obs['celltype_sub'] = adata.obs.index.map(adata_sub.obs['leiden_r0.5'].to_dict())
adata.obs['celltype_sub'].fillna(adata.obs['leiden_r1.5_n2'],inplace=True)

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['celltype_sub'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)

# remove doublet round 3
adata = adata[~(adata.obs['celltype_sub'].isin(['doublet']))].copy()

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['celltype_sub'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)

adata.write('/media/AnalysisFastDisk/NatualCohort_myeloid_Cluster4th_0.7M.h5ad')

genes_keep = [gene for gene in adata.var_names if gene not in remove_genes]
len(genes_keep)
adata_for_deg = adata[:,genes_keep].copy()
rsc.tl.rank_genes_groups_logreg(adata_for_deg, groupby="celltype_sub", use_raw=False)
sc.pl.rank_genes_groups(adata_for_deg, n_genes=20)

# COSG
import cosg as cosg
import importlib
importlib.reload(cosg)

cosg.cosg(adata_for_deg,
    key_added='cosg',
    mu=1,
    n_genes_user=50,
    groupby='celltype_sub')

colnames = ['names', 'scores']
test = [pd.DataFrame(adata_for_deg.uns["cosg"][c]) for c in colnames]
test = pd.concat(test, axis=1, names=[None, 'group'], keys=colnames)

markers = {}
cats = adata_for_deg.obs['celltype_sub'].cat.categories
for i, c in enumerate(cats):
    cell_type_df = test.loc[:, 'names'][str(c)]
    scores_df = test.loc[:, 'scores'][str(c)]
    markers[str(c)] = cell_type_df.values.tolist()[:10]

sc.pl.dotplot(adata_for_deg, 
              var_names = markers, 
              groupby='celltype_sub',
              cmap='Spectral_r',
              standard_scale='var',save='_Myeloid1st_dotplot_celltype_sub.pdf') 

marker_gene=pd.DataFrame(adata_for_deg.uns['cosg']['names'])
marker_gene.to_excel('CIMA/scRNA/Myeloid1st_dotplot_celltype_sub_marker_gene.xlsx')

marker_gene_scores=pd.DataFrame(adata_for_deg.uns['cosg']['scores'])
marker_gene_scores.to_excel('CIMA/scRNA/Myeloid1st_dotplot_celltype_sub_marker_gene_scores.xlsx')

# Annotation
adata.write('CIMA/scRNA/CIMA_myeloid_Annotation1st.h5ad')
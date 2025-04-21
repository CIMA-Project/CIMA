"""
**************************************
*  @Author  :   Yuhui Zheng
*  @Time    :   2024/12/19
*  @Project :   CIMA
*  @FileName:   CIMA_scRNA_Preprocessing.py
*  @Summary :   CIMA scRNA-seq data preprocessing, clustering and macro-cluster annotation。
**************************************
"""

# CIMA scRNA-seq Preprocessing
import scanpy as sc
from multiprocessing import Pool,Lock
import os
import anndata
import numpy as np
import scrublet as scr
import pandas as pd

def read_C4_mtx(sample):
    
    if os.path.isdir(path + sample):

        sub_adata = anndata.read_mtx(path + sample + '/matrix.mtx.gz').T
        features = pd.read_csv(path + sample + '/features.tsv.gz', header=None, sep='\t')
        sub_adata.var_names = anndata.utils.make_index_unique(pd.Index(features[0].values))
        sub_adata.obs_names = pd.read_csv(path + sample + '/barcodes.tsv.gz', header=None)[0].values
        sub_adata.obs_names = sample + '_' + sub_adata.obs_names
        sub_adata.obs['library'] = sample
        sub_adata.obs['sample'] = sample[4:sample.rfind('-')]
        
        scrub = scr.Scrublet(sub_adata.X, expected_doublet_rate=0.06)
        doublet_scores = scrub.scrub_doublets(min_counts=3, min_cells=3,
                                                    log_transform = True,
                                                    min_gene_variability_pctl=85,
                                                    n_prin_comps=30)[0]
        predicted_doublets = scrub.call_doublets(threshold=0.2)

        sub_adata.obs['doublet_scores'] = doublet_scores
        sub_adata.obs['predicted_doublets'] = predicted_doublets
        
        sub_adata.write('CIMA/scRNA/H5adFile/'+sample+'.h5ad')

        
def init(l):
    global lock
    lock = l


path = 'CIMA/scRNA/Matrix/'
samples = os.listdir(path)
lock = Lock()
pool = Pool(initializer=init, initargs=(lock,))
pool.map(read_C4_mtx,samples)
pool.close()
pool.join()



# CIMA scRNA-seq Preprocessing
import scanpy as sc
import anndata
import cupy as cp
import numpy as np
import pandas as pd
import time
import rapids_singlecell as rsc
from rapids_singlecell.cunnData import cunnData
import rapids_scanpy_funcs
import warnings
warnings.filterwarnings("ignore")
from multiprocessing import Pool,Lock
import os

np.random.seed(66)

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=150, facecolor='white')

import rmm
rmm.reinitialize(
    managed_memory=True)
cp.cuda.set_allocator(rmm.rmm_cupy_allocator)

# Multi-threaded reading of h5ad files
def read_h5ad_file(sample):
    sub_adata = anndata.read_h5ad(sample)
    return sub_adata

def init(l):
    global lock
    lock = l

path = 'CIMA/scRNA/H5adFile/'
samples = [path + i for i in os.listdir(path) ]
lock = Lock()
pool = Pool(initializer=init, initargs=(lock,))
sub_adata = pool.map(read_h5ad_file,samples)
pool.close()
pool.join()

adata = anndata.concat(sub_adata, join="outer")

# QC
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['hb'] = adata.var_names.str.contains('^HB[^(P)]')
adata.var['rp'] = adata.var_names.str.match(r'^RP[SL][0-9]')
adata.var['ncRNA'] = adata.var_names.str.match(r'^[A-Z][A-Z][0-9].*\.[0-9]')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','rp','ncRNA','hb'], percent_top=None, log1p = False, inplace=True)
adata.obs['n_counts'] = adata.X.sum(axis=1)

original_num_cells = adata.n_obs
adata = adata[adata.obs['predicted_doublets'] == False, :]
num_filtered_cells = adata.n_obs
num_cells_difference = original_num_cells - num_filtered_cells
print(f"Doublets Filtered out {num_cells_difference} cells. Original number of cells: {original_num_cells}, after filtering: {num_filtered_cells}")
#Doublets Filtered out 296892 cells. Original number of cells: 8031549, after filtering: 7734657

original_num_cells = adata.n_obs
adata = adata[adata.obs['pct_counts_mt'] < 15, :]
num_filtered_cells = adata.n_obs
num_cells_difference = original_num_cells - num_filtered_cells
print(f"MT Filtered out {num_cells_difference} cells. Original number of cells: {original_num_cells}, after filtering: {num_filtered_cells}")
#MT Filtered out 285815 cells. Original number of cells: 7734657, after filtering: 7448842

original_num_cells = adata.n_obs
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_cells(adata, max_genes=6000)
sc.pp.filter_cells(adata, min_counts=1000)
sc.pp.filter_cells(adata, max_counts=25000)
num_filtered_cells = adata.n_obs
num_cells_difference = original_num_cells - num_filtered_cells
print(f"Filtered out {num_cells_difference} cells. Original number of cells: {original_num_cells}, after filtering: {num_filtered_cells}")
#Filtered out 545410 cells. Original number of cells: 7448842, after filtering: 6903432

adata.write('/media/CIMA/CIMA_scRNA_Filtered_6.9M.h5ad')

# HVG selection
adata = sc.read_h5ad('/media/CIMA/CIMA_scRNA_Filtered_6.9M.h5ad')
sc.pp.highly_variable_genes(adata,n_top_genes=2500,flavor="seurat_v3")

mt_genes = list(adata.var_names[adata.var_names.str.match(r'^MT-')])
#hb_genes = list(adata.var_names[adata.var_names.str.contains('^HB[^(P)]')]) 
rp_genes = list(adata.var_names[adata.var_names.str.match(r'^RP[SL]')]) 
ncRNA_genes = list(adata.var_names[adata.var_names.str.match(r'^[A-Z][A-Z][0-9].*\.[0-9]')])
LINC_genes = list(adata.var_names[adata.var_names.str.match(r'(^LOC|LINC)[1-9]*')])
#IG_genes = list(adata.var_names[adata.var_names.str.match(r'^IG[HKL]')]) 
#IG_genes = list(set(IG_genes).difference(set(['IGHM','IGHE','IGHD'])))
#malat1_genes = list(adata.var_names[adata.var_names.str.startswith('MALAT1')])

#remove mt, rp, ncRNA and LINC gene
remove_genes = mt_genes + rp_genes + ncRNA_genes + LINC_genes

hvg_list = adata.var_names[adata.var.highly_variable].to_list()
pd.DataFrame(hvg_list).to_excel('CIMA/scRNA/CIMA_hvg_list.xlsx')
hvg_remove = [gene for gene in hvg_list if gene in remove_genes]
pd.DataFrame(hvg_remove).to_excel('CIMA/scRNA/CIMA_hvg_remove.xlsx')
hvg_keep = [gene for gene in hvg_list if gene not in remove_genes]
pd.DataFrame(hvg_keep).to_excel('CIMA/scRNA/CIMA_hvg_keep.xlsx')

is_keep = np.isin(adata.var_names, hvg_keep)
adata.var['highly_variable'] = is_keep
adata.var.highly_variable.value_counts()

adata.layers["counts"] = adata.X.copy()
adata.write('CIMA/scRNA/CIMA_HVG_6.9M.h5ad')

# Log Nor PCA
adata = sc.read_h5ad('CIMA/scRNA/CIMA_HVG_6.9M.h5ad')
sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)
adata.write('CIMA/scRNA/CIMA_Norlog_6.9M.h5ad')

adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps = 50)
sc.pl.pca_variance_ratio(adata_hvg, log=True, n_pcs=50)
sc.pl.pca_loadings(adata_hvg, components=range(1,30))
adata.uns['pca'] = adata_hvg.uns['pca'].copy()
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()

# Harmony
rsc.tl.harmony_integrate(adata,'sample')
adata.write('CIMA/scRNA/CIMA_hvg_harmony_6.9M.h5ad')

# UMAP
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=19, use_rep='X_pca_harmony')
sc.tl.umap(adata, min_dist=0.3)
sc.tl.leiden(adata, resolution=1.5, n_iterations=2, key_added='leiden_r1.5_n2')
adata.write('CIMA/scRNA/CIMA_UMAP_6.9M.h5ad')
sc.pl.umap(adata,color='leiden_r1.5_n2')

from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
cluster_cols = ["#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", 
                "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", 
                "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D", 
                "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999", 
                "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000", 
                "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"]

sc.pl.umap(adata,color=['HBA1','MZB1','LILRA4','PPBP',
                        'CD3D','CD4','CD8A','NCAM1','CD14','FCGR3A','CD1C','IGHD','CD27','CD34','AIF1','MKI67','CD19','XBP1'])

Annotation_1st = [
    'T cells',#0
    'T cells',
    'T cells',
    'NK cells',
    'T cells',
    'cMono',  #5
    'T cells',
    'T cells',
    'T cells',
    'T cells',
    'Naive B',#10
    
    'NK cells',#11
    'Memory B',
    'NK cells',
    'T cells',
    'cMono',#15
    'T cells',
    'T cells',
    'ncMono',
    'T cells',
    'doublet',#20
    
    'cDC',#21
    'Plasma&Cycling T',
    'pDC',
    'Erythrocyte',
    'Megakaryocyte',
    'HSPC']
adata.obs['Annotation_1st']=adata.obs['leiden_r1.5_n2'].replace(to_replace = list(adata.obs['leiden_r1.5_n2'].cat.categories.values),
                                                         value = Annotation_1st)

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['Annotation_1st'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)
    
adata.write('CIMA/scRNA/CIMA_Annotation_1st_6.9M.h5ad')


#B cell subset
B_cells = adata[adata.obs['celltype_1st'].isin(['Naive B','Memory B','Plasma&Cycling T'])].copy()
B_cells.write('CIMA/scRNA/CIMA_B_cells_0.6M.h5ad')

#Myeloid cell subset
myeloid = adata[adata.obs['celltype_1st'].isin(['pDC','cDC','cMono','ncMono','Megakaryocyte','HSPC'])].copy()
myeloid.write('CIMA/scRNA/CIMA_myeloid_0.9M.h5ad')


#T&NK cell subset
B_Annotation1st_Metadata = pd.read_csv('CIMA/scRNA/B_Annotation1st_Metadata.csv',index_col=0)
adata.obs['celltype_2nd'] = adata.obs['celltype_1st']
adata.obs['celltype_2nd'] = adata.obs.index.map(B_Annotation1st_Metadata['celltype_B_1st'].to_dict())
adata.obs['celltype_2nd'].fillna(adata.obs['celltype_1st'],inplace=True)

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['celltype_2nd'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)
    
myeloid = sc.read_h5ad('CIMA/scRNA/CIMA_myeloid_Cluster3rd_0.7M.h5ad')
adata.obs['celltype_3rd'] = adata.obs['celltype_2nd']
adata.obs['celltype_3rd'] = adata.obs.index.map(myeloid.obs['leiden_r1.5_n2'].to_dict())
adata.obs['celltype_3rd'].fillna(adata.obs['celltype_2nd'],inplace=True)

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, color=['celltype_3rd'], legend_loc='on data',
               legend_fontsize=8, legend_fontoutline=2)

TNK = adata[adata.obs['celltype_3rd'].isin(['NK cells','T cells','Cycling T','19'])].copy()
TNK.write('cima/scRNA/CIMA_TNK_5.2M.h5ad')
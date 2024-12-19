"""
**************************************
*  @Author  :   Yuhui Zheng
*  @Time    :   2024/12/19
*  @Project :   CIMA
*  @FileName:   Metacell_Generation.py
*  @Summary :   Use the SEACells core to assign metacell labels to each cell type in each sample 
                and generate Metacell for scRNA and scATAC.
**************************************
"""

import os
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
import anndata as ad
import numpy as np
import pandas as pd
from SEACells import core
from tqdm import tqdm
from sklearn.preprocessing import normalize
import sklearn.utils.extmath
import scipy.sparse
from typing import Any, Mapping, Optional, TypeVar, Union

Array = Union[np.ndarray, scipy.sparse.spmatrix]

def tfidf(X: Array) -> Array:
    idf = X.shape[0] / X.sum(axis=0)
    if scipy.sparse.issparse(X):
        tf = X.multiply(1 / X.sum(axis=1))
        return tf.multiply(idf)
    else:
        tf = X / X.sum(axis=1, keepdims=True)
        return tf * idf

# Assign metacell labels to each cell type in each sample
def assigning_metacell_labels(
        ad, samples="sample", celltypes='celltype'
):
    sample_list = ad.obs[samples].unique().to_list()
    for sample in tqdm(sample_list, desc=samples):
        temp_sample_adata = ad[ad.obs[samples]==sample]
        
        celltype_list = temp_sample_adata.obs[celltypes].unique().to_list()
        for celltype in tqdm(celltype_list,desc=celltypes):
            temp_celltype_adata = temp_sample_adata[temp_sample_adata.obs[celltypes]==celltype].copy()
        
            if temp_celltype_adata.n_obs < 100:
                temp_celltype_adata.obs['SEACell_ID'] = sample + "-" + celltype + "-SEACell-0"
            else:
                n_SEACells = int(np.floor(temp_celltype_adata.n_obs / 50)) +1
                build_kernel_on = 'X_lsi'
                n_waypoint_eigs = 4
                
                X = tfidf(temp_celltype_adata.X)
                X_norm = normalize(X, norm="l1")
                X_norm = np.log1p(X_norm * 1e4)
                X_lsi = sklearn.utils.extmath.randomized_svd(X_norm, 50, **{"random_state":666})[0]
                X_lsi -= X_lsi.mean(axis=1, keepdims=True)
                X_lsi /= X_lsi.std(axis=1, ddof=1, keepdims=True)
                temp_celltype_adata.obsm["X_lsi"] = X_lsi
                
                model = core.SEACells(temp_celltype_adata, 
                                build_kernel_on=build_kernel_on, 
                                n_SEACells=n_SEACells, 
                                n_waypoint_eigs=n_waypoint_eigs,
                                convergence_epsilon = 1e-5,
                                            use_gpu=False)
                model.construct_kernel_matrix()
                M = model.kernel_matrix
                # Initialize archetypes
                model.initialize_archetypes()
                model.fit(min_iter=10, max_iter=500)
        
                temp_celltype_adata.obs['SEACell_ID'] = sample + "-" + celltype  + '-' + temp_celltype_adata.obs['SEACell']
            
            ad.obs.loc[temp_celltype_adata.obs.index, 'SEACell_ID'] = temp_celltype_adata.obs['SEACell_ID']

    return ad

# Modify the summarize_by_SEACell function of SEACells
def summarize_by_SEACell_modified(
        ad, SEACells_label="SEACell", celltype_label=None, summarize_layer="counts"
):

    import scanpy as sc
    from tqdm import tqdm
    import pandas as pd
    import numpy as np
    from scipy.sparse import csr_matrix

    # Set of metacells
    metacells = ad.obs[SEACells_label].unique()

    # Summary matrix
    summ_matrix = pd.DataFrame(0.0, index=metacells, columns=ad.var_names)

    for m in tqdm(summ_matrix.index):
        cells = ad.obs_names[ad.obs[SEACells_label] == m]
        summ_matrix.loc[m, :] = np.ravel(
            ad[cells, :].X.sum(axis=0)
        )

    # Counts
    meta_ad = sc.AnnData(csr_matrix(summ_matrix), dtype=csr_matrix(summ_matrix).dtype)
    meta_ad.obs_names, meta_ad.var_names = summ_matrix.index.astype(str), ad.var_names

    return meta_ad


# Generate scRNA Metacells
rna = ad.read_h5ad('CIMA/scRNA/CIMA_scRNA_Annotation.h5ad')
rna = assigning_metacell_labels(rna, samples="sample", celltypes='cell_type_l4')

rna.obs['SEACell'] = rna.obs['SEACell_ID'].astype('category')
rna.obs[['SEACell']].to_csv('CIMA/scRNA/CIMA_scRNA_Metacell.csv')

rna_metacell = summarize_by_SEACell_modified(rna)
rna_metacell.write_h5ad('CIMA/scRNA/CIMA_scRNA_Metacell.h5ad')


# Generate scATAC Metacells
atac = ad.read_h5ad('CIMA/scATAC/CIMA_scATAC_Annotation_PeakRecalling_3762242cells_338036peaks.h5ad')
atac = assigning_metacell_labels(atac, samples="sample", celltypes='cell_type_l4')

atac.obs['SEACell'] = atac.obs['SEACell_ID'].astype('category')
atac.obs[['SEACell']].to_csv('CIMA/scATAC/CIMA_scATAC_PeakRecalling_Metacell.csv')

atac_metacell = summarize_by_SEACell_modified(atac)
atac_metacell.write_h5ad('CIMA/scATAC/CIMA_scATAC_PeakRecalling_Metacell.h5ad')

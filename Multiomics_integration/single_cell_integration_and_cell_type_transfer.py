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
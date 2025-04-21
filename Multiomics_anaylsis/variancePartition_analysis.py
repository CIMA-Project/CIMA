#!/usr/bin/env python
# coding: utf-8

# # Generate data

# In[ ]:


# Define a function to process the data for each gene
def process_gene_chunk(chunk, groups=None):
    try:
        # For a given gene block, group and sum by cell type and sample ID
        chunk_sum = chunk.groupby(groups).sum()
        return chunk_sum
    except ValueError:
        print("Please provide grouped keys to process. Try again...")

# pseudobulk_data(DataFrame), each row represents a cell type and sample combination, each column represents a gene, and the value is the sum of the raw counts of all cells in that cell type and sample.
def grouped_obs_sum(adata, group_keys, chunk_size = 1000):
    
    # Initialize an empty DataFrame to store the pseudo data
    pseudobulk_data = pd.DataFrame()

    # Get a list of group combinations
    groups=[adata.obs[key] for key in group_keys]
    
    # Use adata's block processing function
    # The chunk_size here can be adjusted according to your memory size
    chunk_size = chunk_size  # Each time we process n genes
    for start in range(0, adata.shape[1], chunk_size):
        end = min(start + chunk_size, adata.shape[1])
        gene_chunk = adata[:, start:end].to_df()
        
        chunk_sum = process_gene_chunk(gene_chunk, groups)
        # Accumulate the results into pseudobulk_data
        pseudobulk_data = pd.concat([pseudobulk_data, chunk_sum], axis=1)
    return pseudobulk_data


# In[361]:


gene_chunk = adata_rna[:, 1:10].to_df()
groups=[adata_rna.obs[key] for key in ['final_annotation', 'sample']]
pseudobulk_data=process_gene_chunk(gene_chunk, groups)


# In[372]:


pseudobulk_df=grouped_obs_sum(adata_rna, ['sample', 'final_annotation'], chunk_size = 1000)


# In[388]:


pseudobulk_df.index=[i[0] + ':' + i[1] for i in pseudobulk_df.index]
# pseudobulk_df.to_csv('../processed_data/scRNA/pseudobulk_BySampleCelltype.csv')


# In[458]:


pseudobulk_df


# In[442]:


obs=pd.DataFrame([pseudobulk_df.index], index = ['barcode_names'], columns=pseudobulk_df.index).T
obs['individual']=[i.split(':')[0] for i in pseudobulk_df.index]
obs['celltype']=[i.split(':')[1] for i in pseudobulk_df.index]
obs


# In[443]:


var=pd.DataFrame([pseudobulk_df.columns], index = ['gene_name'], columns=pseudobulk_df.columns).T
var


# In[ ]:


###
meta_df=pd.read_excel('../Nature_Cohorts_PBMC_scATAC/20230807_WGS-ATAC-RNA.xlsx')
meta_df=meta_df[['BGE编号', '招募时间', 'RNA实验时间', '录入性别', '录入年龄']]
meta_df.columns=['individual', 'recruitment_time', 'RNA_experiment_time', 'sex', 'age']
meta_df.sex = meta_df['sex'].map({'男':'male', '女':'female'}).astype('category')
meta_df=meta_df.set_index('individual')
meta_df.head()


# In[ ]:


# Convert pseudobulk data back to AnnData format for saving
pseudobulk_adata = sc.AnnData(pseudobulk_df, obs=obs, var=var)

# add additional meta_info
pseudobulk_adata.obs['recruitment_time']=meta_df.loc[pseudobulk_adata.obs['individual']]['recruitment_time'].astype(str).values
pseudobulk_adata.obs['RNA_experiment_time']=meta_df.loc[pseudobulk_adata.obs['individual']]['RNA_experiment_time'].astype(str).values
pseudobulk_adata.obs['age']=meta_df.loc[pseudobulk_adata.obs['individual']]['age'].values
pseudobulk_adata.obs['sex']=meta_df.loc[pseudobulk_adata.obs['individual']]['sex'].values


# In[493]:


pseudobulk_adata.write_h5ad('../processed_data/scRNA/pseudobulk_BySampleCelltype.h5ad', compression="gzip")


# In[494]:


import os
import scipy
from scipy import io
def h5ad_to_mtx(adata, destination=None):
    if not os.path.exists(destination):
        os.makedirs(destination)
    pd.DataFrame(adata.var.index).to_csv(os.path.join(destination, "genes.tsv" ),   sep = "\t", index = False, header=False)
    pd.DataFrame(adata.obs.index).to_csv(os.path.join(destination, "barcodes.tsv"), sep = "\t", index = False, header=False)
    adata.obs.to_csv(os.path.join(destination, "metadata.tsv"), sep = "\t", index = True)
    #io.mmwrite(os.path.join(destination, "matrix.mtx"), scipy.sparse.csr_matrix(adata.X))
    #io.mmwrite(os.path.join(destination, "matrix.mtx"), scipy.sparse.csc_matrix(adata.X))
    io.mmwrite(os.path.join(destination, "matrix.mtx"), scipy.sparse.coo_matrix(adata.X.T))
    return


# In[495]:


h5ad_to_mtx(pseudobulk_adata, destination='../processed_data/scRNA/pseudobulk_data')


# # variancePartition analysis

# In[1]:


library(dreamlet)
library(muscat)
library(ExperimentHub)
library(zenith)
library(scater)
library(Seurat)
library("variancePartition")


# In[6]:


sessionInfo()


# ## pre-processing

# In[18]:


pbmc.data <- Read10X(data.dir = "../processed_data/scRNA/pseudobulk_data/", gene.column = 1)


# In[5]:


meta.data <- read.csv('../processed_data/scRNA/pseudobulk_data/metadata.tsv', sep='\t', row.names = 1)


# In[42]:


pbmc <- CreateSeuratObject(counts = pbmc.data, meta.data = meta.data, min.cells = 0, min.features = 0)


# In[76]:


head(pbmc@meta.data)


# In[81]:


# filter-out celltypes (based on single-cell results)
celltype_df=read.csv('../processed_data/scRNA/ExpressedSamples_MeanCounts_ByCelltype.csv')
celltype_df=subset(celltype_df, sample > 200 & Counts > 9)
#
pbmc <- subset(pbmc, celltype %in% celltype_df$final_annotation)


# In[92]:


# filter-out low-expressed genes
counts <- pbmc@assays$RNA@counts
gene.percent.expressed <- rowMeans(counts > 1)*100
#
genes.filter <- names(gene.percent.expressed[gene.percent.expressed > 10])  #select genes expressed in at least 10% of cells


# In[93]:


length(genes.filter)


# In[97]:


counts.sub <- counts[genes.filter,]
pbmc <- CreateSeuratObject(counts=counts.sub, meta.data = pbmc@meta.data[c('barcode_names', 'individual', 'celltype', 'recruitment_time', 'RNA_experiment_time', 'sex', 'age')],
                                        min.cells = 0, min.features = 0)


# In[101]:


head(pbmc@meta.data)
pbmc


# In[200]:


min(pbmc$nCount_RNA); max(pbmc$nCount_RNA)
min(pbmc$nFeature_RNA); max(pbmc$nFeature_RNA)


# In[206]:


# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0, log = FALSE)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0, log = TRUE)


# In[4]:


saveRDS(pbmc, '../processed_data/scRNA/filtered_pseudobulk_BySampleCelltype.rds')


# * filter-out pseudobulk samples (because some samples don't have cells of subtypes)

# In[4]:


pbmc <- readRDS('../processed_data/scRNA/filtered_pseudobulk_BySampleCelltype.rds')


# In[6]:


pbmc <- subset(pbmc, nFeature_RNA > 5000)


# In[8]:


# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0, log = FALSE)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0, log = TRUE)


# In[9]:


min(pbmc$nCount_RNA); max(pbmc$nCount_RNA)
min(pbmc$nFeature_RNA); max(pbmc$nFeature_RNA)


# ## raw pipe (variancePartition)

# In[4]:


# library(variancePartition)
library(edgeR)
library(BiocParallel)


# In[11]:


# normalize RNA-seq counts
dge <- DGEList(counts = pbmc@assays$RNA@counts)
dge <- calcNormFactors(dge)


# In[32]:


# specify formula with random effect for Individual
#form <- ~  age + (1 | individual) + (1 | celltype) + (1 | sex) + (1 | recruitment_time) + (1 | RNA_experiment_time)
form <- ~  age + (1 | individual) + (1 | celltype) + (1 | sex)

# compute observation weights
#vobj <- voomWithDreamWeights(dge[1:200, ], form, pbmc@meta.data) # test using small gene size
vobj <- voomWithDreamWeights(dge, form, pbmc@meta.data)


# In[34]:


# variancePartition seamlessly deals with the result of voom()
# by default, it seamlessly models the precision weights
# This can be turned off with useWeights=FALSE
varPart <- fitExtractVarPartModel(vobj, form, pbmc@meta.data)


# In[10]:


### V1
#write.csv(varPart, '../processed_data/scRNA/VarPartition_result/varPartResults_15063genes_Dreamlet.csv', quote=FALSE)

### V2
#write.csv(varPart, '../processed_data/scRNA/VarPartition_result/version2/varPartResults_15062genes_Dreamlet.csv', quote=FALSE)
varPart <- read.csv('../processed_data/scRNA/VarPartition_result/version2/varPartResults_15062genes_Dreamlet.csv', row.names = 1)


# In[12]:


vp <- sortCols(varPart, FUN = median)


# In[13]:


plotPercentBars(vp[1:10, ])


# In[14]:


plotVarPart(vp)


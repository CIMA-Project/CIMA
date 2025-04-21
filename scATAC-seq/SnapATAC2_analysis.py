#!/usr/bin/env python
# coding: utf-8

import snapatac2 as snap
snap.__version__

# In[2]:
import sys,os
import pandas as pd
import subprocess
import numpy as np
import scanpy as sc

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

import matplotlib.pyplot as plt
get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'retina'")

# # Identify outliers and doublets (using R)
# In[ ]:
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(patchwork)
library(cowplot)
set.seed(1234)

######################################## load functions ########################################
source('skin_utility.r')
source('seurat_helpers.r')

# In[ ]:
#path1 <- '../Nature_Poputation_project/data/mpy_New'
#path2 <- '../Nature_Poputation_project/data/yy_New'
path3 <- '../Nature_Poputation_project/data/dss_New'

data.list=list.dirs(path3, recursive=FALSE)
files=vapply(strsplit(data.list,"/"), `[`, 9, FUN.VALUE=character(1))
sample_names=gsub('_web_1|_web_0', '', files)

fragment_paths = c()
barcode_paths = c()
peak_paths = c()
metadata_paths = c()
for (i in 1:length(data.list)){
    fragment_paths=c(fragment_paths, paste0(data.list[i], '/', sample_names[i], '.fragments.tsv.gz'))
    barcode_paths =c(barcode_paths, paste0(data.list[i], '/', sample_names[i], '_Peak/barcodes.tsv'))
    peak_paths=c(peak_paths, paste0(data.list[i], '/', sample_names[i], '_Peak/peak.bed'))
    metadata_paths=c(metadata_paths, paste0(data.list[i], '/', sample_names[i], '.Metadata.tsv'))
}

output_path='../Natural_Populations_Blood/QC/qc_SingleLevel/dss_New/'

for (i in 1:length(peak_paths)){
    peaks=Processing_peaks(peak_paths[i])
    peaks=GRangesToString(peaks, sep = c("-", "-"))

    counts <- LoadData(data_path=paste0(data.list[i], '/', sample_names[i], '_Peak'))
    atac.assay <- CreateChromatinAssay(counts, sep = c(":", "-"), fragments = fragment_paths[i])
    seurat.obj <- CreateSeuratObject(atac.assay, assay = "ATAC")

    seurat.obj=seurat.obj[peaks,]
    seurat.obj$orig.ident=sample_names[i]

    seurat.obj <- annotate_ATAC_obj2(seurat.obj, annotation='hg38')

    seurat.obj <- QC_ATAC(seurat.obj, fragment_file = fragment_paths[i], genome = "hg38")

    metadata <- read.csv(file = metadata_paths[i], header = TRUE, row.names = 1, sep='')
    metadata$uniqueCellFrags <- metadata$uniqueFrags + metadata$uniqueMitoFrags
    seurat.obj <- AddMetaData(seurat.obj, metadata = metadata)

    seurat.obj <- runscDblFinder(seurat.obj, assay='ATAC')

    write.csv(seurat.obj@meta.data, paste0(output_path, sample_names[i], '_meta.csv'), quote=FALSE)

    print (paste0('Num of processing files:', i))
}

# ## plot  before QC
# In[ ]:
library(dplyr)
library(ggplot2)
library(patchwork)

is_outlier <- function(dataframe, metric=NULL, nmads=NULL) {
  if (!metric %in% colnames(dataframe)) {
    stop("metric must be a valid column name")
  }
  M <- dataframe[,metric]
  median <- median(M)
  mad <- mad(M)
  outlier <- (M < median - nmads * mad) | (median + nmads * mad < M)
  return(outlier)
}

meta_files1=list.files('../Natural_Populations_Blood/QC/qc_SingleLevel/mpy_New', full.names=TRUE)
meta_files2=list.files('../Natural_Populations_Blood/QC/qc_SingleLevel/yy_New', full.names=TRUE)
meta_files3=list.files('../Natural_Populations_Blood/QC/qc_SingleLevel/dss_New', full.names=TRUE)
meta_files4=list.files('../Natural_Populations_Blood/QC/qc_SingleLevel/old', full.names=TRUE)
meta_files = c(meta_files1, meta_files2, meta_files3, meta_files4)

meta_datas=data.frame()
for (i in 1:length(meta_files)){
    meta_data=read.csv(meta_files[i], row.names=1)
    meta_datas=rbind(meta_datas,meta_data)
}

meta_datas <- subset(meta_datas, nCount_ATAC > 0)

options(repr.plot.width=18, repr.plot.height=10)
p1=ggplot(meta_datas, aes(x=nCount_ATAC)) + geom_density(alpha=.2, fill="#FF6666") + scale_x_continuous(limits = c(0, 20000))
p2=ggplot(meta_datas, aes(x=TSS.enrichment)) + geom_density(alpha=.2, fill="#FF6666") + scale_x_continuous(limits = c(0, 20)) 
p3=ggplot(meta_datas, aes(x=FRiP)) + geom_density(alpha=.2, fill="#FF6666")
p4=ggplot(meta_datas, aes(x=blacklist_fraction)) + geom_density(alpha=.2, fill="#FF6666") + scale_x_continuous(limits = c(0, 0.02)) 
p5=ggplot(meta_datas, aes(x=MitoProportion)) + geom_density(alpha=.2, fill="#FF6666") + scale_x_continuous(limits = c(0, 0.2))
p6=ggplot(meta_datas, aes(x=nucleosome_signal)) + geom_density(alpha=.2, fill="#FF6666") + scale_x_continuous(limits = c(0, 2))
p1+p2+p3+p4+p5+p6

df=meta_datas %>%  group_by(orig.ident) %>%  summarize(mean_Frags = mean(nCount_ATAC), mean_TSS = mean(TSS.enrichment), mean_FRiP = mean(FRiP),
                                                       mean_BF = mean(blacklist_fraction), mean_MP = mean(MitoProportion), mean_NS = mean(nucleosome_signal))
df=df %>% arrange(mean_Frags, mean_FRiP, mean_TSS)

options(repr.plot.width=18, repr.plot.height=8)
par(mfrow=c(2, 3), cex.axis = 1.6, cex.lab = 1.6)
style = theme(axis.text.x = element_text(color = "grey20", size = 24, angle = 90, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 24, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 24, angle = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 24, angle = 90, face = "plain"))
plot(df$mean_Frags, col="red")
plot(df$mean_TSS, col="red")
plot(df$mean_NS, col="red")
plot(df$mean_FRiP, col="green")
plot(df$mean_BF, col="green")
plot(df$mean_MP, col="green")

num_df <- as.data.frame(table(meta_datas$orig.ident))
num_df <- num_df %>% arrange(Freq)
head(num_df)

options(repr.plot.width=7, repr.plot.height=5)
plot(num_df$Freq, col="black")

meta_datas$outlier <- outlier

meta_datas=subset(meta_datas, outlier==FALSE)
meta_datas=subset(meta_datas, scDblFinder.class=='singlet')

options(repr.plot.width=18, repr.plot.height=10)
p1=ggplot(meta_datas, aes(x=nCount_ATAC)) + geom_density(alpha=.2, fill="#FF6666") + scale_x_continuous(limits = c(0, 20000))
p2=ggplot(meta_datas, aes(x=TSS.enrichment)) + geom_density(alpha=.2, fill="#FF6666") + scale_x_continuous(limits = c(0, 20)) 
p3=ggplot(meta_datas, aes(x=FRiP)) + geom_density(alpha=.2, fill="#FF6666")
p4=ggplot(meta_datas, aes(x=blacklist_fraction)) + geom_density(alpha=.2, fill="#FF6666") + scale_x_continuous(limits = c(0, 0.02)) 
p5=ggplot(meta_datas, aes(x=MitoProportion)) + geom_density(alpha=.2, fill="#FF6666") + scale_x_continuous(limits = c(0, 0.2))
p6=ggplot(meta_datas, aes(x=nucleosome_signal)) + geom_density(alpha=.2, fill="#FF6666") + scale_x_continuous(limits = c(0, 2))
p1+p2+p3+p4+p5+p6

df=meta_datas %>%  group_by(orig.ident) %>%  summarize(mean_Frags = mean(nCount_ATAC), mean_TSS = mean(TSS.enrichment), mean_FRiP = mean(FRiP),
                                                       mean_BF = mean(blacklist_fraction), mean_MP = mean(MitoProportion), mean_NS = mean(nucleosome_signal))
df=df %>% arrange(mean_Frags, mean_FRiP, mean_TSS)

write.csv(df, '../Natural_Populations_Blood/QC/statstics_SampleLevel/sample_meta.csv', quote=FALSE)

options(repr.plot.width=18, repr.plot.height=8)
par(mfrow=c(2, 3), cex.axis = 1.6, cex.lab = 1.6)
plot(df$mean_Frags, col="red")
plot(df$mean_TSS, col="red")
plot(df$mean_NS, col="red")
plot(df$mean_FRiP, col="green")
plot(df$mean_BF, col="green")
plot(df$mean_MP, col="green")

num_df <- as.data.frame(table(meta_datas$orig.ident))
num_df <- num_df %>% arrange(Freq)
head(num_df)

options(repr.plot.width=7, repr.plot.height=5)
plot(num_df$Freq, col="black")

meta_df=pd.read_csv('../Nature_Cohorts_PBMC_scATAC/qc_SingleLevel/qc_meta.csv',index_col=0)
meta_df = meta_df[(meta_df['scDblFinder.class']=='singlet') & (meta_df['outlier'] == False)]
meta_df.shape
meta_df = meta_df[(meta_df['nCount_ATAC'] > 1000) & (meta_df['TSS.enrichment'] > 5) & (meta_df['FRiP'] > 0.6) & (meta_df['tssProportion'] > 0.3)]
meta_df.shape

data_dir = "../scATAC/NEW_Data/Fragments/"
output_dir = "../Nature_Cohorts_PBMC_scATAC/h5ad_output"
os.makedirs(output_dir, exist_ok=True)
fragment_files = [f'{data_dir}/{fl}' for fl in os.listdir(data_dir) if fl.endswith(".tsv.gz")]

data_dir = "../Nature_Cohorts_PBMC_scATAC/Modified_Fragments"
fragment_files2 = [f'{data_dir}/{fl}' for fl in os.listdir(data_dir) if fl.endswith(".tsv.gz")]
sample_names2=[f.split('/')[-1] for f in fragment_files2]
fragment_files=[f for f in fragment_files if f.split('/')[-1] not in sample_names2]

fragment_files.extend(fragment_files2)

outputs = []
for fl in fragment_files:
    name = fl.split('/')[-1].split('.tsv.gz')[0]
    outputs.append(f'{output_dir}/{name}.h5ad')

sample_names=[f.split('/')[-1].split('.fragments.tsv.gz')[0] for f in fragment_files]

len(fragment_files), len(outputs), len(sample_names)

cellPass=meta_df.index.tolist()

get_ipython().run_cell_magic('time', '', 'adatas = snap.pp.import_data(\n    fragment_files,\n    file=outputs,\n    chrom_sizes=snap.genome.hg38,\n    min_num_fragments=500,\n    sorted_by_barcode=False,\n    whitelist=cellPass\n)')

snap.pp.add_tile_matrix(adatas, bin_size=5000)
snap.pp.select_features(adatas, n_features=50000)

get_ipython().run_cell_magic('time', '', 'data = snap.AnnDataSet(\n    adatas=[(f.filename.split(\'/\')[-1].split(\'.fragments.h5ad\')[0], f) for f in adatas],\n    filename="./h5ad_output/data.h5ads"\n)')

data.close()

meta_df = pd.read_csv('./Nature_Cohorts_PBMC_scATAC/qc_SingleLevel/qc_meta.csv',index_col=0)
meta_df = meta_df[(meta_df['scDblFinder.class']=='singlet') & (meta_df['outlier'] == False)]
meta_df = meta_df[(meta_df['nCount_ATAC'] > 1000) & (meta_df['TSS.enrichment'] > 5) & (meta_df['FRiP'] > 0.6) & (meta_df['tssProportion'] > 0.3)]

import pandas as pd
from scipy.stats import spearmanr

def calculate_d_statistic(samples):
    n = len(samples)
    d_statistic = pd.Series(index=samples.index)
    
    for i, sample in samples.iterrows():
        correlations = []
        for j, other_sample in samples.iterrows():
            if i != j:
                correlation, _ = spearmanr(sample, other_sample)
                correlations.append(correlation)
        
        d_statistic[i] = pd.Series(correlations).median()
    
    return d_statistic

def identify_outliers(samples, alpha=0.05):
    d_statistic = calculate_d_statistic(samples)
    outlier_threshold = d_statistic.quantile(alpha)
    outliers = d_statistic[d_statistic <= outlier_threshold]
    return outliers

columns=['orig.ident', 'nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal',
         'FRiP', 'blacklist_fraction', 'MitoProportion', 'tssProportion']
df=meta_df[columns].groupby('orig.ident').mean().sort_values('nCount_ATAC')

df['cell_num']=list(meta_df.value_counts(["orig.ident"]).loc[df.index])

outliers1=df[(df['nCount_ATAC'] < 1000) | (df['TSS.enrichment'] < 5) | (df['FRiP'] < 0.6) | (df['tssProportion'] < 0.3)| (df['cell_num'] < 500)]

corr_df=calculate_d_statistic(norm_df)
outliers2=corr_df[corr_df<0.2]

len([i for i in outliers2.index if i in outliers1.index])

remove_samples=['20221115-ATAC-E-B21853979833-1', '20230216-ATAC-E-B21806592250-a-1', '20230216-ATAC-E-B21803652303-a-1', '20230216-ATAC-E-B21803652303-a-2',
                '20230216-ATAC-E-B21806592250-a-2', '20230217-ATAC-E-B21812497385-b-2','20230216-ATAC-E-B21536696547-a-2', '20230209-ATAC-E-B21116522625-b-1',
                '20230217-ATAC-E-B21812497385-b-1', '20230216-ATAC-E-B21536696547-a-1','20230209-ATAC-E-B21116522625-b-2', '20230215-ATAC-E-B21784276925-a-1',
                '20230216-ATAC-E-B21791266764-a-2', '20230110-ATAC-E-B21703926764-a-1','20230302-ATAC-E-B21141858405-b-1', '20230216-ATAC-E-B21830607992-b-2',
                '20230216-ATAC-E-B21792779337-a-1', '20230216-ATAC-E-B21806212813-a-2','20230216-ATAC-E-B21806212813-a-1', '20230216-ATAC-E-B21791266764-a-1',
                '20230215-ATAC-E-B21781969131-a-2', '20230217-ATAC-E-B21182908058-b-2','20230110-ATAC-E-B21703926764-a-2', '20230217-ATAC-E-B21182908058-b-1',
                '20230209-ATAC-E-B21149562361-b-1', '20230215-ATAC-E-B21781969131-a-1','20230215-ATAC-E-B21784276925-a-2', '20230216-ATAC-E-B21792779337-a-2',
                '20221121-ATAC-E-B21446062078-a-2', '20230222-ATAC-E-B21244819307-b-2','20230222-ATAC-E-B21244819307-b-1', '20221121-ATAC-E-B21446062078-a-1',
                '20230221-ATAC-E-B21567528994-a-2', '20230221-ATAC-E-B21229854160-a-1','20230221-ATAC-E-B21567528994-a-1', '20230224-ATAC-E-B21577416048-b-1',
                '20230224-ATAC-E-B21577416048-b-2', '20230221-ATAC-E-B21229962323-a-1','20230110-ATAC-E-B21703620342-a-1', '20230224-ATAC-E-B21574116910-b-2',
                '20230110-ATAC-E-B21700277504-a-1', '20230209-ATAC-E-B21986816169-B-2','20230110-ATAC-E-B21672724169-a-1', '20230110-ATAC-E-B21696276791-a-1',
                '20230110-ATAC-E-B21700277504-a-2', '20230110-ATAC-E-B21696276791-a-2','20230221-ATAC-E-B21564294836-a-2', '20230221-ATAC-E-B21564294836-a-1',
                '20230110-ATAC-E-B21698275649-a-2', '20230209-ATAC-E-B21358717623-b-1','20230110-ATAC-E-B21698275649-a-1', '20230110-ATAC-E-B21681728186-a-2',
                '20230209-ATAC-E-B21358717623-b-2', '20230110-ATAC-E-B21681728186-a-1','20230110-ATAC-E-B21672724169-a-2', '20230221-ATAC-E-B21574133918-a-1',
                '20230110-ATAC-E-B21697788806-a-1', '20230110-ATAC-E-B21697788806-a-2','20230221-ATAC-E-B21574133918-a-2', '20230110-ATAC-E-B21678902021-a-2',
                '20230110-ATAC-E-B21678902021-a-1', '20230110-ATAC-E-B21694489992-a-2','20230110-ATAC-E-B21694489992-a-1']

data_dir='../Nature_Cohorts_PBMC_scATAC/h5ad_output'
h5ad_files = [f'{data_dir}/{fl}' for fl in os.listdir(data_dir) if fl.endswith(".fragments.h5ad")]
h5ad_files=[h for h in h5ad_files if h.split('/')[-1].split('.fragments.h5ad')[0] not in remove_samples]

adatas=[snap.read(i) for i in h5ad_files]

get_ipython().run_cell_magic('time', '', 'data = snap.AnnDataSet(\n    adatas=[(f.filename.split(\'/\')[-1].split(\'.fragments.h5ad\')[0], f) for f in adatas],\n    filename="./h5ad_output/data.h5ads"\n)')

data.obs['library'] = [('-').join(i.split('-')[2:]) for i in data.obs['sample']]
data.obs['batch'] = [i.split('-')[0] for i in data.obs['sample']]
data.obs['individual'] = [('-').join(i.split('-')[2:4]) for i in data.obs['sample']]
data.obs['sample_label'] = [('-').join(i.split('-')[2:]) for i in data.obs['sample']]
data.obs['sample_label'] = [i.replace('-1', '') for i in data.obs['sample_label']]
data.obs['sample_label'] = [i.replace('-2', '') for i in data.obs['sample_label']]

print (len(set(data.obs['sample'])))
print (len(set(data.obs['library'])))
print (len(set(data.obs['batch'])))
print (len(set(data.obs['individual'])))
print (len(set(data.obs['sample_label'])))

print(f'Number of cells: {data.n_obs}')
print(f'Number of unique barcodes: {np.unique(data.obs_names).size}')

get_ipython().run_cell_magic('time', '', 'snap.pp.select_features(data, n_features=200000)\nsnap.tl.spectral(data)')

get_ipython().run_cell_magic('time', '', 'snap.pp.harmony(data, batch="individual", max_iter_harmony=30, key_added=\'X_spectral_harmony\')')

get_ipython().run_cell_magic('time', '', 'snap.tl.umap(data, use_rep="X_spectral_harmony", use_dims=30, random_state=42)')

get_ipython().run_cell_magic('time', '', 'snap.pp.knn(data, use_rep="X_spectral_harmony")')

get_ipython().run_cell_magic('time', '', 'snap.tl.leiden(data, resolution=0.5)')

snap.pl.umap(data, color="leiden", height=600, marker_opacity=1, interactive=False, show=True, scale=1)

snap.pl.umap(data, color="leiden", height=600, marker_size=0.5, interactive=False, show=True, scale=1)

snap.pl.umap(data, color="sample", width=800, height=600, marker_size=0.6, interactive=False, show=True, scale=1)

snap.pl.umap(data, color="individual", width=700, height=600, marker_size=0.6, interactive=False, show=True, scale=1)

snap.pl.umap(data, color="batch", width=700, height=600, marker_size=0.6, interactive=False, show=True, scale=1)

from collections import Counter
Counter(data.obs['leiden'])

cross_df=pd.crosstab(data.obs['leiden'], data.obs['batch'])
cross_df.shape

import seaborn as sns
sns.clustermap(cross_df, standard_scale=0, figsize=(16,8))

cross_df=pd.crosstab(data.obs['leiden'], data.obs['individual'])
cross_df.shape

import seaborn as sns
sns.set(style='whitegrid')
sns.clustermap(cross_df, standard_scale=0)

cross_df=pd.crosstab(data.obs['leiden'], data.obs['sample'])
cross_df=pd.crosstab(data.obs['leiden'], data.obs['library'])

import seaborn as sns
sns.clustermap(cross_df, standard_scale=0)

get_ipython().run_cell_magic('time', '', 'gene_matrix = snap.pp.make_gene_matrix(data, snap.genome.hg38)\ngene_matrix')

gene_matrix.obsm["X_umap"] = data.obsm["X_umap"]

sc.pl.umap(gene_matrix, color=["leiden"], legend_loc='on data')

remove_clusters=['8', '12', '15', '16', '17', '18', '19']
gene_matrix=gene_matrix[~gene_matrix.obs['leiden'].isin(remove_clusters)]

sc.settings.set_figure_params(dpi=80, figsize=(6, 6),facecolor='white')
sc.pl.umap(gene_matrix, color=["leiden"], legend_loc='on data')

gene_matrix.obs['nCount_ATAC']=meta_df.loc[gene_matrix.obs.index]['nCount_ATAC'].copy()
gene_matrix.obs['TSS.enrichment']=meta_df.loc[gene_matrix.obs.index]['TSS.enrichment'].copy()
gene_matrix.obs['FRiP']=meta_df.loc[gene_matrix.obs.index]['FRiP'].copy()

sc.settings.set_figure_params(dpi=80, figsize=(6, 6),facecolor='white')
sc.pl.umap(gene_matrix, color=["nCount_ATAC",'TSS.enrichment', 'FRiP'])

gene_matrix.write("./h5ad_output/pbmc_gene_mat.h5ad", compression='gzip')

df=pd.read_excel('./20230807_WGS-ATAC-RNA.xlsx')
df=df[~df['final processed'].isin(['RNA and ATAC are both removed', 'SNP, RNA and ATAC were eliminated'])]

overlap_samples=[i for i in set(gene_matrix.obs['individual']) if i in set(df['BGEnumber'])]
len(overlap_samples)

cellID = gene_matrix.obs[gene_matrix.obs['individual'].isin(overlap_samples)].index.tolist()
len(cellID)

gene_matrix2=gene_matrix[cellID, ]

gene_matrix2.write("./SnapATAC2/pbmc_filtered_genescore.h5ad", compression="gzip")

sub_data=data.subset(obs_indices=cellID, var_indices=None, out='./SnapATAC2/final_filtered', backend=None)

sub_data[0].close()

data.close()

meta_data=pd.read_csv('../Natural_Cohort/scATAC/SnapATAC2/final_filtered/meta_scATAC_3762242cells.csv', index_col=0)

data_dir = "../scATAC/NEW_Data/Fragments"
output_dir = "../Nature_Cohorts_PBMC_scATAC/h5ad_output"
os.makedirs(output_dir, exist_ok=True)
fragment_files = [f'{data_dir}/{fl}' for fl in os.listdir(data_dir) if fl.endswith(".tsv.gz")]

data_dir = "../Nature_Cohorts_PBMC_scATAC/Modified_Fragments"
fragment_files2 = [f'{data_dir}/{fl}' for fl in os.listdir(data_dir) if fl.endswith(".tsv.gz")]

sample_names2=[f.split('/')[-1] for f in fragment_files2]
fragment_files=[f for f in fragment_files if f.split('/')[-1] not in sample_names2]

fragment_files.extend(fragment_files2)

sample_names=meta_data['sample'].unique()
fragment_files=[f for f in fragment_files if f.split('/')[-1].split('.fragments.tsv.gz')[0] in sample_names]

outputs = []
for fl in fragment_files:
    name = fl.split('/')[-1].split('.tsv.gz')[0]
    outputs.append(f'{output_dir}/{name}.h5ad')

cellPass = meta_data.index.tolist()

get_ipython().run_cell_magic('time', '', 'adatas = snap.pp.import_data(\n    fragment_files,\n    file=outputs,\n    chrom_sizes=snap.genome.hg38,\n    min_num_fragments=500,\n    sorted_by_barcode=False,\n    whitelist=cellPass\n)')

data = snap.read_dataset("./h5ad_output/data.h5ads")

adata = data.to_adata()

import sys

def get_size_in_units(obj, unit='MB'):
    size_in_bytes = sys.getsizeof(obj)
    
    if unit == 'KB':
        size = size_in_bytes / 1024
    elif unit == 'MB':
        size = size_in_bytes / (1024 ** 2)
    elif unit == 'GB':
        size = size_in_bytes / (1024 ** 3)
    else:
        raise ValueError("Unit must be 'KB', 'MB', or 'GB'")
    
    return size
print(f"Memory size of list_var: {get_size_in_units(adata, 'GB')} GB")

adata.write_h5ad('./h5ad_output/adata_5120bp.h5ad', compression='gzip')

data = snap.read_dataset("./h5ad_output/data.h5ads")

adata = data.to_adata()

print(f"Memory size of list_var: {get_size_in_units(adata, 'GB')} GB")

adata.write_h5ad('./h5ad_output/adata_10240bp.h5ad', compression='gzip')

data = snap.read_dataset("./h5ad_output/data.h5ads")

adata = data.to_adata()

print(f"Memory size of list_var: {get_size_in_units(adata, 'GB')} GB")

adata.write_h5ad('./h5ad_output/adata_51200bp.h5ad', compression='gzip')


# In[ ]:





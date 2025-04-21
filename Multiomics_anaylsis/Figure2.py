#!/usr/bin/env python
# coding: utf-8

# # Fig. 2B-PCA-Py
# In[1]:
import scanpy as sc
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, figsize=(4, 4),facecolor='white')
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
import matplotlib.pyplot as pl
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import seaborn as sns
sns.set_style(style='ticks')
sc.set_figure_params(dpi=80,dpi_save=300,vector_friendly="pdf")
# In[4]:
#metadta
meta_df=pd.read_excel('../bulk_data/Table/CIMA_Table_S1.xlsx', index_col=0)
meta_df
# In[9]:
#cell ratio
ratio_df=pd.read_excel('../bulk_data/Table/CIMA_Table_S1.xlsx', index_col=0, sheet_name=6)
ratio_df=ratio_df.T
ratio_df
# In[10]:
def norm_zero2one(df):
    normalized_df=(df-df.min())/(df.max()-df.min())
    return normalized_df
# In[11]:
#PCA
ratio_df=norm_zero2one(ratio_df)
adata = sc.AnnData(X=ratio_df, obs=meta_df.loc[ratio_df.index])
sc.pl.highest_expr_genes(adata, n_top=10, )
sc.pp.highly_variable_genes(adata, min_mean=0.1, max_mean=1, min_disp=0.5)
adata.var[adata.var['highly_variable']==True]
sc.pl.highly_variable_genes(adata)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)
sc.pl.pca(adata, color=['Age','sex'],size=200, annotate_var_explained=True, frameon=True, ncols=2,
          #save="RNA_proportion_PC.pdf"
         )
# # Fig. 2D-Boxplot-R
# In[ ]:
library(ggplot2)
library(future)
library(tidyverse)
library(pheatmap)
library(ggpubr)
library(ggchicklet)
library(magrittr)
library(ggh4x)
library(rstatix)
library(ggsignif)
library(ggnewscale)
library(patchwork)
library(gapminder)
library(ggprism)
# In[ ]:
#data
tmp_data <- read.csv('../1_Natural_Cohort/bulkData/Data/metabolome.csv',check.names = FALSE,header=T)
tmp_data <- data.frame(t(tmp_data))
colnames(tmp_data) <- tmp_data[1,]
tmp_data <- tmp_data[-1,]
tmp_data
# In[ ]:
#metadata
tmp_meatdata <- read.csv('../1_Natural_Cohort/bulkData/Data/sample_level_meta.csv',check.names = FALSE,header=T,row.names = 1)
tmp_meatdata <- tmp_meatdata[,c(1,5)]
tmp_meatdata <- na.omit(tmp_meatdata)
tmp_meatdata
# In[ ]:
#merge
common_cols <- intersect(rownames(tmp_data), rownames(tmp_meatdata))
tmp_data_sub <- tmp_data[common_cols,]
tmp_meatdata_sub <- tmp_meatdata[common_cols,]
tmp_data_sub
# In[ ]:
#correlation
mytheme <- theme_prism(base_family="",base_fontface="plain") + 
  theme(strip.text = element_text(size = 10,angle=0,vjust = 0.5,hjust = 0.5),
        axis.line = element_line(color = "black",size = 0.2),
        axis.ticks = element_line(size = 0.2),
        axis.text.y = element_text(color = "black",size = 12),
        axis.text.x = element_text(color = "black",size = 12),
        axis.title = element_text(color = "black",size = 12),
        legend.position = "none")

data <- cbind(tmp_data_sub,tmp_meatdata_sub)
rownames(data) <- rownames(tmp_meatdata_sub)
for (i in 1:321){
    for (j in 323){
        subdata <- data[,c(i,j)]
        subdata[, 1] <- as.numeric(as.character(subdata[, 1]))
        subdata[, 1] <- (subdata[, 1] - min(subdata[, 1])) / (max(subdata[, 1]) - min(subdata[, 1]))
        p <- ggplot(subdata, aes(subdata[,2],subdata[,1]))+RotatedAxis()+
        stat_boxplot(aes(color=sex),geom ="errorbar",width=0.2,linewidth=0.25)+
        geom_boxplot(aes(color=sex),outlier.shape = NA,linewidth=0.25)+
        geom_point(aes(color=sex),size=1.5,alpha=0.5,position = "jitter")+
        scale_color_manual(values = c("#BB2F29","#142F67"))+
        stat_compare_means(comparisons = list(c("Female","Male")),
                     method="wilcox.test",
                     paired = F,
                     size=2.5,
                     symnum.args = list(cutpoints = c(0,0.001,0.01,0.05,1),
                                        symbols = c("***","**","*","ns")))+
        labs(title = "")+
        theme(plot.title = element_text(size=8,hjust=0.5))+
        xlab("")+ylab(paste0(colnames(subdata)[1]))+ mytheme
        ggsave(paste0("../1_Natural_Cohort/bulkData/figures/boxplot/metabolome/",colnames(subdata)[1],"_",colnames(subdata)[2]," boxplot.pdf"),p,width=2.1,height=2.8)
        }
}

# # Fig. 2E-ScatterPlot-R
# In[ ]:
library(ggplot2)
library(future)
library(tidyverse)
library(pheatmap)
library(ggpubr)
library(ggchicklet)
library(magrittr)
library(ggh4x)
library(rstatix)
library(ggsignif)
library(ggnewscale)
library(patchwork)
library(gapminder)
library(ggprism)
# In[ ]:
#data
tmp_data <- read.csv('../1_Natural_Cohort/bulkData/Data/celltype_proportion_RNA.csv',check.names = FALSE,header=T)
tmp_data <- data.frame(t(tmp_data))
colnames(tmp_data) <- tmp_data[1,]
tmp_data <- tmp_data[-1,]
tmp_data
# In[ ]:
#metadata
tmp_meatdata <- read.csv('../1_Natural_Cohort/bulkData/Data/sample_level_meta.csv',check.names = FALSE,header=T,row.names = 1)
tmp_meatdata <- tmp_meatdata[,c(1,5)]
tmp_meatdata$sex <- ifelse(tmp_meatdata$sex == 'Female', 0, 1)
tmp_meatdata <- na.omit(tmp_meatdata)
tmp_meatdata
# In[ ]:
#merge
common_cols <- intersect(rownames(tmp_data), rownames(tmp_meatdata))
tmp_data_sub <- tmp_data[common_cols,]
tmp_meatdata_sub <- tmp_meatdata[common_cols,]
# In[ ]:
#correlation
mytheme <- theme_prism(base_family="",base_fontface="plain") + 
  theme(strip.text = element_text(size = 10,angle=0,vjust = 0.5,hjust = 0.5),
        axis.line = element_line(color = "black",size = 0.2),
        axis.ticks = element_line(size = 0.2),
        axis.text.y = element_text(color = "black",size = 12),
        axis.text.x = element_text(color = "black",size = 12),
        axis.title = element_text(color = "black",size = 12),
        legend.position = "none")

data <- cbind(tmp_data_sub,tmp_meatdata_sub)
rownames(data) <- rownames(tmp_meatdata_sub)

for (i in 1:73){
    for (j in 75){
        subdata <- data[,c(i,j)]
        subdata[, 1] <- as.numeric(as.character(subdata[, 1]))
        celltype <- colnames(subdata)[1]
        colnames(subdata)[1] <- "Proportion"
        p <- ggplot(data=subdata, aes(x=subdata[,2],y=subdata[,1]))+
        geom_point(aes(color=Age),size=1.5,alpha=0.75,position = "jitter")+
        scale_color_gradient(low = "#A7C8DF", high = "#142F67") +
        stat_smooth(method="lm",se=T)+
        stat_cor(data=subdata, method = "spearman")+
        labs(x=paste0(colnames(subdata)[2]),y=paste0(colnames(subdata)[1]))+ggtitle(paste0(celltype," spearman correlation"))+ mytheme
        ggsave(paste0("../1_Natural_Cohort/bulkData/figures/scatterplot/ATAC proportion/",celltype,"_",colnames(subdata)[2]," correlation.pdf"),p,width=4,height=4)
        }
}


# # Fig. 2F-Heatmap-R

# In[ ]:
library(ggplot2)
library(future)
library(tidyverse)
library(pheatmap)
# In[ ]:
#metadata
metadata <- as.data.frame(t(read.csv("../1_Natural_Cohort/bulkData/Data/sample_level_meta.csv",row.names = 1,header = T,check.names = F)))
metadata <- as.data.frame(t(metadata))
metadata <- metadata[,c(1,4,5)]
metadata
# In[ ]:
#MOFA_model_factors
data1 <- as.data.frame(t(read.csv("../1_Natural_Cohort/MOFA/Data/RNA_ATAC_2group_model_factors.csv",row.names = 1,header = T,check.names = F)))
data1 <- as.data.frame(t(data1))
data1
# In[ ]:
#merge
common_cols <- intersect(rownames(metadata), rownames(data1))
metadata <- metadata[common_cols,]
data1 <- data1[common_cols,]
tmp <- cbind(metadata,data1)
tmp$sex <- ifelse(tmp$sex == "Female", 0, 1)
tmp[,1] <- as.numeric(as.character(tmp[,1]))
tmp[,2] <- as.numeric(as.character(tmp[,2]))
tmp[,3] <- as.numeric(as.character(tmp[,3]))
tmp
# In[ ]:
#correlation
correlation_results <- data.frame()
for (i in 1:3) {
  for (j in 4:length(colnames(tmp))) {
    cor_test <- cor.test(tmp[[i]], tmp[[j]])
    correlation_results <- rbind(correlation_results, data.frame(
      Group = names(tmp)[i],
      Factor = names(tmp)[j],
      Correlation = cor_test$estimate,
      P_value = cor_test$p.value
    ))
  }
}

#change data format
Correlation <- correlation_results[,c(1,2,3)]
P_value <- correlation_results[,c(1,2,4)]
Correlation_new <- split(Correlation, Correlation$Group)
Correlation_new <- data.frame(lapply(Correlation_new, function(x) x$Correlation))
rownames(Correlation_new) <- unique(Correlation$Factor)
P_value_new <- split(P_value, P_value$Group)
P_value_new <- data.frame(lapply(P_value_new, function(x) x$P_value))
rownames(P_value_new) <- unique(P_value$Factor)
P_value_new <- -log(P_value_new,10)
P_value_new <- t(P_value_new)
P_value_new
                                 
#pheatmap
p <- pheatmap(P_value_new,
              cluster_cols = F,
              cluster_rows = F,
              fontsize=7,
              color = colorRampPalette(c("white","#DEA1A5"))(100), 
              fontsize_col = 8,
              fontsize_row = 8,
              show_colnames = T,
              cellwidth = 12, 
              cellheight = 12,
              # Set title:
              main = "sample_para_2group_model_factors",
              # Set whether the label is displayed
              annotation_legend	= T,
              scale="none",
              border= F)
ggsave("../1_Natural_Cohort/MOFA/Figure/sample_para_2group_model_factors.pdf",p,width=5,height=3,
      limitsize = FALSE)


# # Fig. 2G
# ## Fig. 2G-Heatmap-R
# In[ ]:
library(ggplot2)
library(future)
library(tidyverse)
library(pheatmap)
# In[ ]:
#data
data <- as.data.frame(read.csv("../1_Natural_Cohort/MOFA/Data/RNA_ATAC_2group_model_r2.csv",header = T,check.names = F))
data
# In[ ]:
#2batch
data_list <- split(data, data$Group)
data_group1 <- as.data.frame(data_list[[1]])[,-1]
data_group2 <- as.data.frame(data_list[[2]])[,-1]
# In[ ]:
#top10
subdata_RNA_count <- subdata_RNA_count[c(1:10),]
# In[ ]:
#batch1
data_group1_list <- split(data_group1, data_group1$Factor)
data_group1_new <- data.frame(lapply(data_group1_list, function(x) x$R2))
rownames(data_group1_new) <- unique(data_group1$View)
Factor7 <- data_group1_new[as.character(subdata_RNA_count$Group),c(7)]
data_group1_new <- data.frame(rownames=c(subdata_RNA_count$Group),Factor7=Factor7)
rownames(data_group1_new) <- data_group1_new$rownames
data_group1_new <- subset(data_group1_new, select = -rownames)
data_group1_new
# In[ ]:
#batch2
data_group2_list <- split(data_group2, data_group2$Factor)
data_group2_new <- data.frame(lapply(data_group2_list, function(x) x$R2))
rownames(data_group2_new) <- unique(data_group2$View)
Factor7 <- data_group2_new[as.character(subdata_RNA_count$Group),c(7)]
data_group2_new <- data.frame(rownames=c(subdata_RNA_count$Group),Factor7=Factor7)
rownames(data_group2_new) <- data_group2_new$rownames
data_group2_new <- subset(data_group2_new, select = -rownames)
data_group2_new
# In[ ]:
p <- pheatmap(data_group1_new,
              cluster_cols = F,
              cluster_rows = F,
              fontsize=7,
              color = colorRampPalette(c("white","black"))(100), 
              fontsize_col = 8,
              fontsize_row = 8,
              show_colnames = T,
              cellwidth = 12, 
              cellheight = 12,
              main = "RNA_ATAC_2group_model_r2_group1_top20",
              annotation_legend	= T,
              scale="none",
              border= F)
ggsave("../1_Natural_Cohort/MOFA/Figure/RNA_ATAC_2group_model_r2_group1_top20.pdf",p,width=5,height=18,limitsize = FALSE)


# In[ ]:
p <- pheatmap(data_group2_new,
              cluster_cols = F,
              cluster_rows = F,
              fontsize=7,
              color = colorRampPalette(c("white","black"))(100), 
              fontsize_col = 8,
              fontsize_row = 8,
              show_colnames = T,
              cellwidth = 12, 
              cellheight = 12,
              main = "RNA_ATAC_2group_model_r2_group2_top20",
              annotation_legend	= T,
              scale="none",
              border= F)
ggsave("../1_Natural_Cohort/MOFA/Figure/RNA_ATAC_2group_model_r2_group2_top20.pdf",p,width=5,height=18,limitsize = FALSE)
# ## Fig. 2G-lollipop-R
# In[ ]:
library(ggplot2)
library(future)
library(tidyverse)
library(pheatmap)
# In[ ]:
#data
data <- as.data.frame(read.csv("../1_Natural_Cohort/MOFA/Data/RNA_ATAC_2group_model_weights.csv",header = T,check.names = F))
data <- data[,c(1,8)]
data <- data[order(data$Factor7,decreasing = TRUE), ]
colnames(data) <- c("Name","Factor7")
data$Group <- data$Name
data$Group <- gsub("^[^_]*_", "", data$Group)
data$omic <- sub(".*_", "", data$Group)
data
# In[ ]:
subdata  <- subset(data, abs(Factor7)>3)
subdata
# In[ ]:
#RNA-top10
subdata_RNA <- subset(subdata,omic %in% c("RNA"))
subdata_RNA_count <- data.frame(table(subdata_RNA$Group))
colnames(subdata_RNA_count) <- c("Group","subdata_count")
subdata_RNA_count <- subdata_RNA_count[order(subdata_RNA_count$subdata_count,decreasing = TRUE), ]
subdata_RNA_count$Group <- factor(subdata_RNA_count$Group,levels = subdata_RNA_count$Group)
subdata_RNA_count <- subdata_RNA_count[c(1:10),]
subdata_RNA_count
# In[ ]:
#lollipop
mytheme <- theme_prism(base_family="",base_fontface="plain") + 
  theme(strip.text = element_text(size = 8,angle=0,vjust = 0.5,hjust = 0.5),
        axis.line = element_line(color = "black",size = 0.2),
        axis.ticks = element_line(size = 0.2),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10,angle=30,vjust = 1,hjust = 1),
        axis.title = element_text(color = "black",size = 10),
        legend.position = "none")

p <- ggplot(subdata_RNA_count, aes(x = subdata_count, y = Group)) +
            geom_point(size = 5) + coord_cartesian(xlim = c(40, NA)) +
            geom_bar(stat = "identity", width = 0.1) +
            ggtitle("Number of features RNA") + mytheme
ggsave("../1_Natural_Cohort/MOFA/Figure/RNA_ATAC_2group_model_weights_top.pdf",p,width=5,height=5,limitsize = FALSE)


# # Fig. 2H-Heatmap-R

# In[ ]:
library(Scillus)
library(tidyverse)
library(Seurat)
library(magrittr)
# In[ ]:
pbmc <- readRDS('../processed_data/scRNA/filtered_pseudobulk_BySampleCelltype.rds')
# filter-out pseudobulk samples with expressed genes less than n
pbmc <- subset(pbmc, nFeature_RNA > 5000)
# convert int to string
pbmc@meta.data$recruitment_time <- as.character(pbmc@meta.data[, 'recruitment_time'])
pbmc@meta.data$RNA_experiment_time <- as.character(pbmc@meta.data[, 'RNA_experiment_time'])
# In[ ]:
pos_genes = c('CX3CR1','CCR2','CCR5','GIMAP7','P2RY13','GBP1','CXorf21',
 'HSPA1A','CD180','PRF1','NFE2','FCGR3A','SLC40A1','TLR10',
 'C21orf91','CD200R1','CCR1','GBP4','CCR6','PIGM','PHF23')

neg_genes = c('MAP3K7CL', 'ATP5MGL', 'AL121944.1',
 'NRARP', 'LINC01619', 'AL118516.1', 'AL136454.1',
 'MAFF', 'ATF3', 'IL1B', 'EREG', 'PF4', 'AC245014.3',
 'CCL3', 'G0S2', 'PDE4D', 'DDIT3', 'AREG', 'CXCL8',
 'AL355075.4', 'SLC11A1', 'AC253572.2', 'RGCC',
 'NR4A3', 'PFKFB3', 'IRS2', 'IER3','SERTAD1')


# In[ ]:
pbmc_sub <- subset(pbmc, celltype %in% c('cDC2_CD1C'))
pbmc_sub <- NormalizeData(pbmc_sub, scale.factor = 1e4)
pbmc_sub <- ScaleData(pbmc_sub)
# In[ ]:
avg_exp <- AverageExpression(pbmc_sub, group.by='age', return.seurat = TRUE)
# In[ ]:
avg_exp@meta.data$age=as.integer(row.names(avg_exp@meta.data))
# In[ ]:
continuous_var <- avg_exp$age
custom_intervals <- c(10, 30, 40, 50, 60, 70, 80)
# Use the cut() function to convert continuous variables into categorical variables.
grouped_var <- cut(continuous_var, breaks = custom_intervals, labels = c("20_29", "30_39", "40_49", "50_59", "60_69", "70_79"))
avg_exp$age_group <- grouped_var
# In[ ]:
avg_exp@meta.data
# In[ ]:
options(repr.plot.width=8, repr.plot.height=4)
#pdf(file = "../Figures/Figure2/cDC2_CD1C_PosGenes-AgeTrend.pdf", width = 8, height = 4)
plot_heatmap(dataset = avg_exp, row_font_size = 10,
              markers = pos_genes,
              sort_var = c("age", "age_group"),
              anno_var = c("age", "age_group"),
              anno_colors = list("Reds", # RColorBrewer palette
                                 "Set1"),
             hm_limit = c(-2, 0, 2), #hm_colors = c("blue","white","red")
             )
#dev.off()
# In[ ]:
options(repr.plot.width=8, repr.plot.height=5)
#pdf(file = "../Figures/Figure2/cDC2_CD1C_NegGenes-AgeTrend_ascending.pdf", width = 8, height = 5)
plot_heatmap(dataset = avg_exp, row_font_size = 9,
              markers = neg_genes,
              sort_var = c("age", "age_group"),
              anno_var = c("age", "age_group"),
              anno_colors = list("Reds", # RColorBrewer palette
                                 "Set2"),
             hm_limit = c(-2, 0, 2), #hm_colors = c("purple","black","yellow")
             )
#dev.off()
# # Fig. 2I-2J-2K
# In[ ]:
varPart = pd.read_csv('../processed_data/scRNA/VarPartition_result/varPartResults_15063genes_Dreamlet.csv', index_col=0)
# In[ ]:
# GO_BP_2023
go_dict_BP={}
with open('../downloaded/GeneSets/GO_Biological_Process_2023.txt') as f:
    data=f.readlines()
    data=[i.strip().split('\t') for i in data]
for term in data:
    go_dict_BP[term[0]]=term[2:]
    
# immune-related genesets
import json
with open('../downloaded/GeneSets/20240425_immue_related_geneset_all_368.json') as json_data:
    go_dict_immune = json.load(json_data)

# Hallmark geneset (MsigDB)
go_path='../Natural_Cohort/downloaded/GeneSets/Hallmark_geneset_for_GSVA_analysis'
go_dict_hallmark={}
for term in os.listdir(go_path):
    tmp_df=pd.read_csv(os.path.join(go_path, term))
    term_name=tmp_df.columns[0]
    go_dict_hallmark[term_name] = tmp_df[term_name].to_list()[1:]

# other specific genesets
go_path = '../Natural_Cohort/downloaded/GeneSets/specific_geneset'
go_dict={}

# Aging
SASP_path=os.path.join(go_path, 'SenMayo_SASP.xlsx')
SASP_genes=pd.read_excel(SASP_path)['Gene(human)'].values # human
go_dict['SenMayo_SASP']=SASP_genes

SASP_path=os.path.join(go_path, 'seneset_Tao_CellMetabolism.txt')
go_dict['seneset']=pd.read_csv(SASP_path, header=None)[0].values

# Housekeeping
geneset=pd.read_csv(os.path.join(go_path, 'HOUNKPE_HOUSEKEEPING_GENES.v2023.2.Hs.tsv'), sep='\t', index_col=0).loc['GENE_SYMBOLS', 'HOUNKPE_HOUSEKEEPING_GENES']
go_dict['housekeeping']=geneset.split(',')

# Cell-cycle
gene_path=os.path.join(go_path, 'regev_lab_cell_cycle_genes.txt')
cell_cycle_genes = [x.strip() for x in open(gene_path)]
go_dict['cellcycle']=cell_cycle_genes

geneset=pd.read_csv(os.path.join(go_path, 'BIOCARTA_CELLCYCLE_PATHWAY.v2023.2.Hs.tsv'), sep='\t', index_col=0).loc['GENE_SYMBOLS', 'BIOCARTA_CELLCYCLE_PATHWAY']
go_dict['cellcycle_BIOCARTA']=geneset.split(',')

# Inflammation
#geneset=pd.read_csv(os.path.join(go_path, 'GOBP_INFLAMMATORY_RESPONSE.v2023.2.Hs.tsv'), sep='\t', index_col=0).loc['GENE_SYMBOLS', 'GOBP_INFLAMMATORY_RESPONSE']
#go_dict['inflammation_GOBP']=geneset.split(',')
#
geneset=pd.read_csv(os.path.join(go_path, 'HALLMARK_INFLAMMATORY_RESPONSE.txt'), header=None)[0]
go_dict['inflammation_HALLMARK']=geneset.values

# Drug respnse
geneset=pd.read_csv(os.path.join(go_path, 'GOBP_RESPONSE_TO_DRUG.v2023.2.Hs.tsv'), sep='\t', index_col=0).loc['GENE_SYMBOLS', 'GOBP_RESPONSE_TO_DRUG']
go_dict['DrugResponse_GOBP']=geneset.split(',')

# cytokine and receptors
geneset=pd.read_csv(os.path.join(go_path, 'Cytokines.txt'), sep='\t', index_col=0)['Symbol']
go_dict['Cytokines']=geneset.values

geneset=pd.read_csv(os.path.join(go_path, 'Cytokine_Receptors.txt'), sep='\t', index_col=0)['Symbol']
go_dict['Cytokine_Receptors']=geneset.values

# Chemokines and receptors
geneset=pd.read_csv(os.path.join(go_path, 'Chemokines.txt'), sep='\t', index_col=0)['Symbol']
go_dict['Chemokines']=geneset.values

geneset=pd.read_csv(os.path.join(go_path, 'Chemokine_Receptors.txt'), sep='\t', index_col=0)['Symbol']
go_dict['Chemokine_Receptors']=geneset.values

# Sex difference
df=pd.read_csv(os.path.join(go_path, 'Male_female_diff_pantissue.csv'), index_col=0)
geneset=[i.split('_')[0] for i in df[df['Whole_Blood']!=0].index.to_list()]
go_dict['Sex_diff']=geneset

# Clock gene
geneset=pd.read_csv(os.path.join(go_path, 'core_clock_genes.txt'), sep='\t', header=None)[0].tolist()
go_dict['clock_genes']=geneset


# ## 2I-Barplots

# In[ ]:
varPart_genesets={}
#
for term in go_dict_BP.keys():
    genes=[i for i in go_dict_BP[term] if i in varPart.index.values]
    print (len(genes))
    if len(genes) > 0: ##Select terms with more than 0 genes
        varPart_genesets[term]=varPart.loc[genes].mean()        
#
for term in go_dict_immune.keys():
    genes=[i for i in go_dict_immune[term] if i in varPart.index.values]
    print (len(genes))
    if len(genes) > 0: 
        varPart_genesets[term]=varPart.loc[genes].mean()        
#
for term in go_dict_hallmark.keys():
    genes=[i for i in go_dict_hallmark[term] if i in varPart.index.values]
    print (len(genes))
    if len(genes) > 0: 
        varPart_genesets[term]=varPart.loc[genes].mean()
#
for term in go_dict.keys():
    genes=[i for i in go_dict[term] if i in varPart.index.values]
    print (len(genes))
    if len(genes) > 0: 
        varPart_genesets[term]=varPart.loc[genes].mean()
#
varPart_genesets=pd.DataFrame(varPart_genesets).T
varPart_genesets=varPart_genesets.dropna()
varPart_genesets.head()
#
class_genesets=[]
for i in varPart_genesets.index.values:
    if i in go_dict_BP.keys():
        class_genesets.append('BP')
    elif i in go_dict_immune.keys():
        class_genesets.append('immune')
    elif i in go_dict_hallmark.keys():
        class_genesets.append('hallmark')
    else:
        class_genesets.append('spec')
varPart_genesets['class']=class_genesets
#


# In[ ]:


df_sorted = varPart_genesets.sort_values('individual', ascending=False).reset_index()
# Initialize an empty list to store the rows that meet the conditions
filtered_rows = []

# Record the value of the previous Residuals
previous_residual = 0.3
# Traverse the sorted DataFrame
for index, row in df_sorted.iterrows():
    current_residual = row['celltype']
    # If the current row is greater than or equal to the previous row, keep the row
    if current_residual > previous_residual - 0.01 and current_residual < previous_residual + 0.02:
        #print (True)
        filtered_rows.append(row)
        previous_residual = current_residual

# Create a new DataFrame to store rows that meet the criteria
filtered_df = pd.DataFrame(filtered_rows)
#
filtered_df = pd.concat([df_sorted[0:5], filtered_df]).sort_values('individual', ascending=False).drop_duplicates()
# In[ ]:
ax=filtered_df[['index', 'individual']].plot(x='index', kind='bar', stacked=True, width=1, linewidth=0.1,
                 title='Percent variation of different resources', grid=False, figsize=(4, 3), legend=None)
ax.set_xticks([])
plt.ylabel('Percent of variation of individual', fontsize=12)
#plt.legend(bbox_to_anchor=(1, 0), loc=3, borderaxespad=0)
#plt.savefig('../Figures/Figure2/varPart/varPart-individual_GOterms_barplot.pdf')
plt.show()
# In[ ]:
ax=filtered_df[['index', 'celltype']].plot(x='index', kind='bar', stacked=True, width=1, linewidth=0.1, colormap='Set2',
                 title='Percent variation of different resources', grid=False, figsize=(4, 3), legend=None)
ax.set_xticks([])
plt.ylabel('Percent of variation of celltype', fontsize=12)
#plt.legend(bbox_to_anchor=(1, 0), loc=3, borderaxespad=0)
#plt.savefig('../Figures/Figure2/varPart/varPart-celltype_GOterms_barplot.pdf')
plt.show()
# ## 2J-Boxplots
# In[ ]:
go_dict_varPart={}
merged_df=pd.DataFrame()
#
for term in go_dict.keys():
    genes=[i for i in go_dict[term] if i in varPart.index.values]
    tmp_df=varPart.loc[genes]
    tmp_df['GeneSets']=term
    go_dict_varPart[term]=tmp_df
    #
    merged_df=pd.concat([merged_df, tmp_df])

# all genes
tmp_df=varPart.copy()
tmp_df['GeneSets']='all genes'
go_dict_varPart['all genes']=tmp_df
merged_df=pd.concat([merged_df, tmp_df])


# In[ ]:
median_merged_df=pd.DataFrame(merged_df.groupby('GeneSets')['celltype'].median()).sort_values('celltype')
median_merged_df
# In[ ]:
#
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
# Use Seaborn's boxplot to plot variables
plt.figure(figsize=(5, 5))
sns.boxplot(merged_df, x='GeneSets', y='celltype', order=median_merged_df.index, color='lightcoral', fliersize=0.5)
plt.xticks(rotation=90, fontsize=12)
# Setting the title and axis labels
plt.title('Transcriptome Variation of functional genesets', fontsize=14)
plt.xlabel('genesets', fontsize=12)
plt.ylabel('Percent of variation of celltype', fontsize=12)
plt.axhline(y=0.416482, color='k', linestyle='--', linewidth=0.8)
# Adjust layout
plt.tight_layout()
# save
plt.savefig('../Figures/Figure2/varPart/VarPart-celltype_SpecTerms_boxplot.pdf')
# show
plt.show()
# In[ ]:
median_merged_df=pd.DataFrame(merged_df.groupby('GeneSets')['individual'].median()).sort_values('individual')
median_merged_df
# In[ ]:
# Use Seaborn's boxplot to plot variables
plt.figure(figsize=(5, 5))
sns.boxplot(merged_df, x='GeneSets', y='individual', order=median_merged_df.index, color='lightblue', fliersize=1)
plt.xticks(rotation=90, fontsize=12)
#plt.ylim([0, 1])
# Setting the title and axis labels
plt.title('Transcriptome Variation of functional genesets', fontsize=14)
plt.xlabel('genesets', fontsize=12)
plt.ylabel('Percent of variation of individual', fontsize=12)
plt.axhline(y=0.067043, color='k', linestyle='--', linewidth=0.8)
# Adjust layout
plt.tight_layout()
plt.savefig('../Figures/Figure2/varPart/VarPart-individual_SpecTerms_boxplot.pdf')
plt.show()
# ## 2K-Barplots
# In[ ]:
gene_list=['RPS4Y1', 'XIST', 'PPP1R15A', 'PER1', 'IER5', 'PTPRC', 'UBB', 'CX3CR1', 'IKZF1', 'CCR2', 
           'GIMAP7', 'CCR5', 'BCL2', 'CD180', 'CD33', 'FLT3']
df=varPart.loc[gene_list].reset_index()
df=df.sort_values('individual', ascending=True)
df
# In[ ]:
# Create two subplots, arranged symmetrically in the middle
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, df.shape[0]*0.25))

# The first bar graph shows scATAC data using a logarithmic scale.
ax1.barh(df['index'], df['individual'], color='lightblue')
ax1.set_title('Percent variation of individual', fontsize=10)
ax1.set_xlabel('Percent variation', fontsize=10)
ax1.set_xlim(0, 1)
#ax1.set_xticks([0.1, 0.2, 0.3, 0.4, 0.8, 1])
#ax1.set_xticklabels([0.1, 0.2, 0.3, 0.4, 0.8, 1])
ax1.yaxis.tick_right() # Move the y-axis scale to the right
ax1.invert_xaxis()  # Reverse the x-axis so the bars run from right to left
ax1.tick_params(axis='both', which='major', labelsize=10, labelrotation=0)
ax1.tick_params(axis='both', which='minor', labelsize=10, labelrotation=0)

# The second bar chart shows scRNA data
ax2.barh(df['index'], df['celltype'], color='lightcoral')
ax2.set_title('Percent variation of celltype', fontsize=10)
ax2.set_xlabel('Percent variation', fontsize=10)
ax2.set_xlim(0, 1)
ax2.set_yticklabels([]) # Set the y-axis label to empty, since it is already displayed on the first child
ax2.tick_params(axis='both', which='major', labelsize=10, labelrotation=0)
ax2.tick_params(axis='both', which='minor', labelsize=10, labelrotation=0)
# Adjust sub-image spacing
#plt.subplots_adjust(wspace=0.35)
plt.tight_layout()
plt.savefig('../Figures/Figure2/varPart/varPart_genes-interest_barplot.pdf')
plt.show()
# ## 2K-UMAP plots
# In[ ]:
# annotation of scRNA
adata_rna = sc.read_h5ad('../scRNA_Data/NatualCohort_All_Annotation_Final_reUMAP.h5ad', backed="r")
# In[ ]:
# loading clean data (selecting partion of samples)
clean_df=pd.read_csv('../Natural_Cohort/meta_data/clean_age.csv')
samples=clean_df['Questionnaire'].tolist()
# In[ ]:
sc.settings.set_figure_params(dpi=80, figsize=(4, 4),facecolor='white')
sc.pl.umap(adata_rna[adata_rna.obs['sample'].isin(samples), ], color=['PER1', 'CX3CR1', 'PPP1R15A'], # color_map='magma',
           cmap='Reds', save='UMAP_genes.pdf')


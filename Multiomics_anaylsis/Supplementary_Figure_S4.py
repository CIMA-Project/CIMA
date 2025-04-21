#!/usr/bin/env python
# coding: utf-8

# # 3A-B (R-based)

# In[ ]:
library(Seurat)
library(ggplot2)
library(future)
library(tidyverse)
library(harmony)
library(SeuratDisk)
library(pheatmap)
library(ggpubr)
library(ggchicklet)
library(ggsci)
library(magrittr)
library(ggh4x)
library(rstatix)
library(ggsignif)
library(ggnewscale)
library(patchwork)
library(gapminder)
library(ggprism)

# data <- as.data.frame(read.csv("../1_Natural_Cohort/MOFA/Data/sample_para_2group_model_weights.csv",header = T,check.names = F))
# data <- data[,c(1,3)]
# data <- data[order(data$Factor2,decreasing = TRUE), ]
# colnames(data) <- c("Name","Factor")
# data$Group <- data$Name
# data$Group <- sub(".*_", "", data$Group)
# data

# In[ ]:
colnames(data) <- c("Name","Factor2")
data$Group <- data$Name
data$Group <- sub(".*_", "", data$Group)
data

# In[ ]:
# 5%
data$Factor_abs <- abs(data$Factor2)
data <- data[order(data$Factor_abs,decreasing = TRUE), ]
sig <- round(length(rownames(data))*0.05)
subdata <-  data[c(1:sig),]
subdata

# In[ ]:
subdata_count <- data.frame(table(subdata$Group))
colnames(subdata_count) <- c("Group","subdata_count")
subdata_count

# In[ ]:
data_count <- data.frame(table(data$Group))
colnames(data_count) <- c("Group","data_count")
data_count

# In[ ]:
count <- merge(subdata_count,data_count)
count$percentage <- count$subdata_count/count$data_count
count$Group <- factor(count$Group,levels=c("metabolome","lipidome","RNA","ATAC","biochemical"))
count

# In[ ]:
mytheme <- theme_prism(base_family="",base_fontface="plain") + 
  theme(strip.text = element_text(size = 8,angle=0,vjust = 0.5,hjust = 0.5),
        axis.line = element_line(color = "black",size = 0.2),
        axis.ticks = element_line(size = 0.2),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 10),
        legend.position = "none")
p1 <- ggplot(count, aes(x = subdata_count, y = Group)) +
            geom_bar(stat = "identity", width = 1) +
            ggtitle("Number of features")+mytheme
p2 <- ggplot(count, aes(x = percentage, y = Group)) +
            geom_bar(stat = "identity", width = 1) +
            ggtitle("Top-ranked features(%)")+mytheme
p <- p1|p2
p
ggsave("../1_Natural_Cohort/MOFA/Figure/sample_para_2group_model_weights_top5.pdf",p,width=10,height=5,limitsize = FALSE)

# # 3D (R-based)

# In[ ]:
library(ggplot2)
library(future)
library(tidyverse)
library(pheatmap)

# In[ ]:
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
subdata_RNA <- subset(subdata,omic %in% c("RNA"))
subdata_RNA

# In[ ]:
subdata_RNA_count <- data.frame(table(subdata_RNA$Group))
colnames(subdata_RNA_count) <- c("Group","subdata_count")
subdata_RNA_count <- subdata_RNA_count[order(subdata_RNA_count$subdata_count,decreasing = TRUE), ]
subdata_RNA_count$Group <- factor(subdata_RNA_count$Group,levels = subdata_RNA_count$Group)
subdata_RNA_count

# In[ ]:
mytheme <- theme_prism(base_family="",base_fontface="plain") + 
  theme(strip.text = element_text(size = 8,angle=0,vjust = 0.5,hjust = 0.5),
        axis.line = element_line(color = "black",size = 0.2),
        axis.ticks = element_line(size = 0.2),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10,angle=30,vjust = 1,hjust = 1),
        axis.title = element_text(color = "black",size = 10),
        legend.position = "none")
p1 <- ggplot(subdata_RNA_count, aes(x = Group, y = subdata_count)) +
            coord_cartesian(ylim = c(20, NA)) +
            geom_bar(stat = "identity", width = 0.5) +
            ggtitle("Number of features RNA") + mytheme
p1
ggsave("../1_Natural_Cohort/MOFA/Figure/RNA_ATAC_2group_model_weights_top.pdf",p1,width=12,height=4,limitsize = FALSE)

# # 3E

# In[3]:
import pandas as pd
import numpy as np

# In[1]:
scRNA_l3_annotations={
    'B':['Bn_TCL1A',
        'Switched_Bm_IGHDneg',
        'Transitional_B_SOX4',
        'Unswitched_Bm_CD1C',
        'Switched_Bm_IGHE',
        'Switched_activated_Bm_CD86',
        'pre-Switched_Bm_JAM3',
        'Atypical_Bm_ITGAX',
        'Bn_IL6',
        'Bn_IFIT3',
        'Plasma_IGHA1',
        'Unswitched_Bm_IL6',
        'Plasma_IGHG1',
        'Plasmablast_MKI67',
        'pre-Switched_Bm_IFIT3'
    ],
    'CD4_T':[
        'CD4_Tn_CCR7',
        'CD4_Tem_CCR7neg',
        'CD4_Tfh-like_CXCR5',
        'CD4_Th1-like_GZMK',
        'CD4_Tcm_CXCR5',
        'CD4_CTL_GZMH',
        'CD4_Th17-like_RORC',
        'CD4_Treg_FCRL3',
        'CD4_Treg_FOXP3',
        'CD4_Th_LMNA',
        'CD4_Th_TNFRSF11A',
        'CD4_Tn_SOX4',
        'CD4_Tcm_IFI44L',
        'CD4_Tn_CXCR5',
        'CD4_Th22-like_CCR10',
        'CD4_Th_CCR4',
        'CD4_Tn_LIMS1',
        'CD4_Tr1-like_IL10',
        'CD4_Th_CR1',
        'CD4_Tem_CCR5',
    ],
    'CD8_T':[
        'CD8_Tn_CCR7',
        'CD8_CTL_GZMB',
        'CD8_Tem_CCR7neg',
        'CD8_CTL_GZMK',
        'CD8_Tcm_IFI44L',
        'CD8_Tn_SOX4',
        'CD8_CTL_IFI44L',
    ],
    'Myeloid':[
        'cMono_CD14',
        'ncMono_FCGR3A',
        'cMono_IFI44L',
        'cDC2_CD1C',
        'ncMono_C1QA',
        'cMono_IL1B',
        'pDC_IRF4',
        'intMono_GFRA2',
        'MK_GP9',
        'ncMono_IFIT1',
        'cMono_CXCL10',
        'cDC_CSF2RA',
        'HSPC_CD34',
        'cDC1_BATF3',
        'AS_DC'
    ],
    'NK&ILC':[
        'Mature_NK_dim_FCGR3A',
        'Terminal_NK_dim_CD160neg',
        'Transitional_NK_GZMK',
        'NK_bright_XCL1',
        'Cycling_NK_MKI67',
        'ILC2_IL2RA',
        'Inflamed_NK_dim_IFIT1'
    ],
    'unconvensional_T':[
        'NKT_NCR1',
        'MAIT_SLC4A10',
        'gdT2_GZMH',
        'gdT2_IL12RB2',
        'gdT2_GZMK',
        'Cycling_T_MKI67',
        'gdT1_TRDV1',
        'pre-T-like_CABP4',
        'NKT_IFNG'
    ]
}

# In[ ]:
# loading data
mofa_df=pd.read_excel('../Natural_Cohort/processed_data/MOFA/RNA_ATAC_2group_model_weights.xlsx', index_col=0)
mofa_df['feature']=[i.split('_')[0] for i in mofa_df.index]
mofa_df['omics']=[i.split('_')[-1] for i in mofa_df.index]
mofa_df['celltype']=[('_').join(i.split('_')[1:-1]) for i in mofa_df.index]
mofa_df

## positive-correlated with factor-7

# In[28]:
df=mofa_df[(mofa_df['Factor7'] > 3) & (mofa_df['omics'] == 'RNA')]

# In[29]:
celltype_l1_annotation=[]
for i in df['celltype'].tolist():
    for j in scRNA_l3_annotations.keys():
        if i in scRNA_l3_annotations[j]:
            celltype_l1_annotation.append(j)
df.loc[:,'l1_annotation']=celltype_l1_annotation

# In[31]:
freq_df=df['feature'].value_counts().reset_index()

# In[32]:
celltype_num=[]
for gene in freq_df['feature'].tolist():
    count_=df[df['feature']==gene]['l1_annotation'].value_counts()
    celltype_num.append(len(count_.index.tolist()))
freq_df['celltype_num']=celltype_num

# In[33]:
freq_df['count_modify']=[i if i < 10 else 10 for i in freq_df['count']]

# In[34]:
count_df=freq_df['count_modify'].value_counts().reset_index(name='Frequency')
count_df=count_df.sort_values('count_modify')
count_df

# In[ ]:
# plot PieChart of ovlp and non-ovlp
import matplotlib.pyplot as plt

labels = count_df['count_modify'].values
sizes = count_df['Frequency'].values

colors = plt.cm.tab20.colors
plt.figure(figsize=(6, 6)) 
plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90, textprops={'fontsize': 10})
plt.axis('equal')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=8)
plt.title('Counts of Gene Frequency', fontsize=14)
plt.savefig('../Figures/Supple_Figure2/MOFA/PieChart_Gene-Frequency_counts_MOFA.pdf')
plt.show()

## negative-correlated with factor-7

# In[21]:
df=mofa_df[(mofa_df['Factor7'] < -3) & (mofa_df['omics'] == 'RNA')]

# In[22]:
freq_df=df['feature'].value_counts().reset_index()

# In[23]:
freq_df['count_modify']=[i if i < 10 else 10 for i in freq_df['count']]

# In[24]:
freq_df.to_csv('../processed_data/MOFA/self-processed/NegGene-Factor7_frequency.csv')

# In[25]:
count_df=freq_df['count_modify'].value_counts().reset_index(name='Frequency')
count_df=count_df.sort_values('count_modify')
count_df

# In[ ]:
# plot PieChart of ovlp and non-ovlp
import matplotlib.pyplot as plt

labels = count_df['count_modify'].values
sizes = count_df['Frequency'].values

colors = plt.cm.tab20.colors

plt.figure(figsize=(6, 6))
plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90, textprops={'fontsize': 10})
plt.axis('equal')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=8)

plt.title('Counts of Gene Frequency', fontsize=14)
plt.show()

# # 3F

# In[ ]:
import gseapy as gp
from gseapy import barplot, dotplot

# In[ ]:
plot_df1=pd.read_csv('../processed_data/MOFA/self-processed/PosGene-Factor7_frequency.csv', index_col=0)
plot_df2=pd.read_csv('../processed_data/MOFA/self-processed/NegGene-Factor7_frequency.csv', index_col=0)

# In[ ]:
spec_gene1=plot_df1[plot_df1['count_modify']==1]['feature'].tolist()
shared_gene1=plot_df1[plot_df1['count_modify'] > 4]['feature'].tolist()
spec_gene2=plot_df2[plot_df2['count_modify']==1]['feature'].tolist()
shared_gene2=plot_df2[plot_df2['count_modify'] > 4]['feature'].tolist()

# In[ ]:
len(spec_gene1), len(shared_gene1), len(spec_gene2), len(shared_gene2)

# In[ ]:
df=pd.DataFrame([spec_gene1, shared_gene1, spec_gene2, shared_gene2], index=['pos-specific', 'pos-shared', 'neg-specific', 'neg-shared']).T
df.to_csv('../processed_data/MOFA/self-processed/Genelist_4classes.csv', index=None)

# In[ ]:
enr1 = gp.enrichr(gene_list=spec_gene1, gene_sets=['MSigDB_Hallmark_2020','KEGG_2021_Human','GO_Biological_Process_2021'], organism='human', background=None, outdir=None)
enr2 = gp.enrichr(gene_list=shared_gene1, gene_sets=['MSigDB_Hallmark_2020','KEGG_2021_Human','GO_Biological_Process_2021'], organism='human', background=None, outdir=None)
enr3 = gp.enrichr(gene_list=spec_gene2, gene_sets=['MSigDB_Hallmark_2020','KEGG_2021_Human','GO_Biological_Process_2021'], organism='human', background=None, outdir=None)
enr4 = gp.enrichr(gene_list=shared_gene2, gene_sets=['MSigDB_Hallmark_2020','KEGG_2021_Human','GO_Biological_Process_2021'], organism='human', background=None, outdir=None)

# In[ ]:
enr1.results['obj']='pos-specific'
enr2.results['obj']='pos-shared'
enr3.results['obj']='neg-specific'
enr4.results['obj']='neg-shared'
results=pd.concat([enr1.results, enr2.results, enr3.results, enr4.results])
results1=results[results['Gene_set']=='MSigDB_Hallmark_2020']
results2=results[results['Gene_set']=='KEGG_2021_Human']
results3=results[results['Gene_set']=='GO_Biological_Process_2021']

# In[ ]:
results.to_csv('../processed_data/MOFA/self-processed/GO_enrichment_4classes.csv')
results=pd.read_csv('../processed_data/MOFA/self-processed/GO_enrichment_4classes.csv', index_col=0)

# In[581]:
terms=['positive regulation of histone ubiquitination (GO:0033184)',
       'positive regulation of vascular associated smooth muscle cell migration (GO:1904754)',
       'positive regulation of membrane potential (GO:0045838)',
       'dendritic cell chemotaxis (GO:0002407)',
       'regulation of macrophage migration (GO:1905521)',
       'regulation of T cell mediated immunity (GO:0002709)',
       'protein targeting to ER (GO:0045047)',
       'Ribosome',
       'regulation of p38MAPK cascade (GO:1900744)',
       'regulation of cell adhesion molecule production (GO:0060353)'
       ]

# In[582]:
plot_df=pd.read_csv('../processed_data/MOFA/self-processed/GO_enrichment_4classes.csv', index_col=0)
plot_df=plot_df[plot_df['Term'].isin(terms)]

# In[412]:
ax = dotplot(plot_df,
              column="Adjusted P-value",
              x='obj', 
              cutoff=0.05,
              size=12,
              top_term=5,
              figsize=(2.5, 4), cmap='plasma_r',
              title = "Function enrichment",
              xticklabels_rot=45, 
              show_ring=True, 
              marker='o',
              ofname='../Figures/Supple_Figure2/MOFA/GO_enrichment_4classes.pdf'
             )

# # 3G (R-based)

# In[7]:
library(Scillus)
library(tidyverse)
library(Seurat)
library(magrittr)

# In[36]:
pbmc <- readRDS('../processed_data/scRNA/filtered_pseudobulk_BySampleCelltype.rds')
pbmc <- subset(pbmc, nFeature_RNA > 5000)
pbmc@meta.data$recruitment_time <- as.character(pbmc@meta.data[, 'recruitment_time'])
pbmc@meta.data$RNA_experiment_time <- as.character(pbmc@meta.data[, 'RNA_experiment_time'])

# In[8]:
pbmc <- NormalizeData(pbmc, scale.factor = 1e6)

# In[13]:
pbmc <- ScaleData(pbmc)

# In[126]:
celltypes=c('CD4_Tem_CCR7neg', 'MAIT_SLC4A10', 'CD8_CTL_GZMB', 'cDC2_CD1C', 'MK_GP9',
           'CD4_Tn_CCR7', 'ncMono_FCGR3A', 'CD4_Th17-like_RORC', 'Bn_TCL1A', 'gdT2_IL12RB2')
pos_shared_genes = c('CX3CR1', 'CISH', 'CETN3', 'FCGR3A', 'C14orf119', 'CCR2', 'PRF1')
neg_shared_genes = c('SERTAD1','NR4A2','IER3','MAP3K7CL','PER1','GNG11','AREG','S100A9','S100A8','ATP5MGL','GRASP')

# In[127]:
pbmc_sub <- subset(pbmc, celltype %in% celltypes)
pbmc_sub <- NormalizeData(pbmc_sub, scale.factor = 1e4)
pbmc_sub <- ScaleData(pbmc_sub)

# In[142]:
avg_exp <- AverageExpression(pbmc_sub, group.by=c('celltype', 'age'), return.seurat = TRUE)

# In[163]:
avg_exp@meta.data$age=as.integer(gsub("^.*_", "", row.names(avg_exp@meta.data)))
avg_exp@meta.data$celltype=gsub('[_][^_]+$', '', row.names(avg_exp@meta.data))

# In[ ]:
continuous_var <- avg_exp$age
custom_intervals <- c(10, 30, 40, 50, 60, 70, 80)
grouped_var <- cut(continuous_var, breaks = custom_intervals, labels = c("20_29", "30_39", "40_49", "50_59", "60_69", "70_79"))
avg_exp$age_group <- grouped_var

# In[191]:
for (i in 1:length(celltypes)) {
    pbmc_sub <- subset(pbmc, celltype %in% celltypes[i])
    pbmc_sub <- NormalizeData(pbmc_sub, scale.factor = 1e4)
    pbmc_sub <- ScaleData(pbmc_sub)
    avg_exp <- AverageExpression(pbmc_sub, group.by='age', return.seurat = TRUE)
    avg_exp@meta.data$age=as.integer(row.names(avg_exp@meta.data))
    continuous_var <- avg_exp$age
    custom_intervals <- c(10, 30, 40, 50, 60, 70, 80)
    grouped_var <- cut(continuous_var, breaks = custom_intervals, labels = c("20_29", "30_39", "40_49", "50_59", "60_69", "70_79"))
    avg_exp$age_group <- grouped_var
    
    pdf(file = paste0("../Figures/Figure2/", celltypes[i], "_PosSharedGenes-AgeTrend.pdf"), width = 8, height = 2.5)
    p=plot_heatmap(dataset = avg_exp, row_font_size = 10,
              markers = pos_shared_genes,
              sort_var = c("age", "age_group"),
              anno_var = c("age", "age_group"),
              anno_colors = list("Reds", "Set1"),
              hm_limit = c(-2, 0, 2)
             )
    dev.off()
}

# In[192]:
for (i in 1:length(celltypes)) {
    pbmc_sub <- subset(pbmc, celltype %in% celltypes[i])
    pbmc_sub <- NormalizeData(pbmc_sub, scale.factor = 1e4)
    pbmc_sub <- ScaleData(pbmc_sub)
    avg_exp <- AverageExpression(pbmc_sub, group.by='age', return.seurat = TRUE)
    avg_exp@meta.data$age=as.integer(row.names(avg_exp@meta.data))
    continuous_var <- avg_exp$age
    custom_intervals <- c(10, 30, 40, 50, 60, 70, 80)
    grouped_var <- cut(continuous_var, breaks = custom_intervals, labels = c("20_29", "30_39", "40_49", "50_59", "60_69", "70_79"))
    avg_exp$age_group <- grouped_var
    
    pdf(file = paste0("../Figures/Figure2/MOFA/", celltypes[i], "_NegSharedGenes-AgeTrend.pdf"), width = 8, height = 2.5)
    p=plot_heatmap(dataset = avg_exp, row_font_size = 10,
              markers = neg_shared_genes,
              sort_var = c("age", "age_group"),
              anno_var = c("age", "age_group"),
              anno_colors = list("Reds", "Set1"),
              hm_limit = c(-2, 0, 2)
             )
    dev.off()
}
```
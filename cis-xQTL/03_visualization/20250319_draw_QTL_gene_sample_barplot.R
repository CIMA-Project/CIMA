library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(glue)

plotting_sort <- read.csv('/CIMA/Data/69_celltype_sort.csv',row.names = 1)
L1_L4_mapping <- read.csv('/CIMA/Data/L1_L4.csv')
eQTL_df = read.csv('/CIMA/Data/20241230_xQTL_cell_sample_number/20241230_celltype_70_scRNA.csv')
eQTL_df$sample_1_count = 0
eQTL_df$gene_count = 0

for (celltype in eQTL_df$final_annotation){
  phenotype_df <- fread(glue('/CIMA/Data/eQTL/pseudobulk/{celltype}.csv'))
  eQTL_df[eQTL_df$final_annotation == celltype,"sample_1_count"] <- nrow(phenotype_df)
  eQTL_df[eQTL_df$final_annotation == celltype,"gene_count"] <- (ncol(phenotype_df)-1)
  rm(phenotype_df)
}

eQTL_df <-merge(eQTL_df,L1_L4_mapping,by.x = 'final_annotation',by.y = 'celltype' )
eQTL_df$final_annotation <- factor(eQTL_df$final_annotation,levels = plotting_sort$celltype)
eQTL_df$category =  factor(eQTL_df$category , levels = c("CD4 T", "CD8 T & unconvensional T", "NK&ILC", "B", "Myeloid",'HSPC'))

'''
NK #cca9cf
B #f8b945
Myeloid #91bd4a
HSPC #3ab37b
CD4T #237b9f
CD8T #8ed1e5
'''

# 创建柱形图_gene_counts
barplot <- ggplot(eQTL_df, aes(
  x = final_annotation,             
  y = gene_count,  
  fill = category)) +
  geom_bar(stat = 'identity')+                                # 画柱形图                                             # x轴与y轴互换位置
  geom_text(                                                  # 在图形上加上数字标签
    aes(label=gene_count, hjust = 1 ,angle = 90),
    size=2                                                    # 标签大小
  ) + theme_classic()+
  labs(x = "celltype", y = "Number of Genes in QTL detection")+
  scale_fill_manual(values = c('#237b9f','#8ed1e5','#cca9cf','#f8b945','#91bd4a','#3ab37b'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "None")
barplot

ggsave("/CIMA/Result/plot/basic_number/Genes_eQTL_detection.pdf", plot = barplot, width = 12, height = 4, dpi = 300, device = "pdf")

# 创建柱形图_sample_counts
barplot <- ggplot(eQTL_df, aes(
  x = final_annotation,             
  y = sample_1_count,  
  fill = category)) +
  geom_bar(stat = 'identity')+                                # 画柱形图                                             # x轴与y轴互换位置
  geom_text(                                                  # 在图形上加上数字标签
    aes(label=sample_1_count, hjust = 1 ,angle = 90),
    size=2                                                    # 标签大小
  ) + theme_classic()+
  labs(x = "celltype", y = "Number of Genes in QTL detection")+
  scale_fill_manual(values = c('#237b9f','#8ed1e5','#cca9cf','#f8b945','#91bd4a','#3ab37b'))+
  theme(axis.text.x  = element_blank(),legend.position = "None")
barplot

ggsave("/CIMA/Result/plot/basic_number/samples_eQTL_detection.pdf", plot = barplot, width = 12, height = 1, dpi = 300, device = "pdf")



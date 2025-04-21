library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(glue)

eQTL_all_df = read.csv('/CIMA/Result/summary/20250110_trans_eQTL_fdr005.csv')
CT_full_list = data.frame("V1" = unique(eQTL_all_df$celltype))
CT_full_list$eGene_num <- 0

eGenelist <- list()
for (celltype in CT_full_list$V1){
  print(celltype)
  eGene_df <- eQTL_all_df[eQTL_all_df$celltype == celltype,]
  eGene_number <- length(unique(eGene_df$phenotype_id))
  if (eGene_number > 0) {eGenelist[[celltype]] <- c(unique(eGene_df$phenotype_id))}
  CT_full_list[CT_full_list$V1 == celltype,]$eGene_num <- eGene_number
}

all_eGene <- unique(eQTL_all_df$phenotype_id)

count_gene_occurrences <- function(gene) {
  # 使用sapply检查每个基因是否在R列表的每个元素中出现，并返回逻辑向量
  gene_occurrences <- sapply(eGenelist, function(x) gene %in% x)
  # 使用sum计算逻辑向量中为TRUE的元素个数，即基因在列表中出现的次数
  count <- sum(gene_occurrences)
  #names(count) <- gene
  return(count)
}

# 使用lapply对基因列表中的每个基因应用count_gene_occurrences函数
all_eGene_cscount <- lapply(all_eGene, count_gene_occurrences)
all_eGene_cscount_df = as.data.frame.vector(all_eGene_cscount,row.names = all_eGene)
all_eGene_cscount_df$all_eGene_cscount <- as.numeric(all_eGene_cscount_df$all_eGene_cscount)
table(all_eGene_cscount_df$all_eGene_cscount == 1)

#各个细胞类型中的cs_eGene list的数量
cs_eGenelist <- list()
for (celltype in names(eGenelist)){
  eGenelist_cal <- eGenelist
  eGene_cs <- eGenelist_cal[[celltype]] 
  eGenelist_cal[celltype] <- NULL
  all_rest_eGene <- Reduce(union, eGenelist_cal)
  eGene_cs <- setdiff(eGene_cs, all_rest_eGene)
  cs_eGenelist[[celltype]] <- eGene_cs
}

CT_full_list$cs_eGene_num <- 0
for (celltype in names(cs_eGenelist)){
  CT_full_list[CT_full_list$V1 == celltype,"cs_eGene_num"] <- length(cs_eGenelist[[celltype]])
}


write.csv(CT_full_list,file = '/CIMA/Result/summary/20250323_trans_eQTL_number_by_celltype.csv')

#画饼图
all_eGene_cscount_df = as.data.frame.vector(all_eGene_cscount,row.names = all_eGene)
all_eGene_cscount_df$all_eGene_cscount <- as.numeric(all_eGene_cscount_df$all_eGene_cscount)          
all_eGene_cscount_df_pie <- all_eGene_cscount_df

library(dplyr)

all_eGene_cscount_df_pie <- all_eGene_cscount_df_pie %>%
  mutate(all_eGene_cscount = case_when(
    all_eGene_cscount > 10 ~ "over_10",
    all_eGene_cscount > 1 & all_eGene_cscount <= 10 ~ "2-10",
    all_eGene_cscount == 1 ~ "1",
    TRUE ~ as.character(all_eGene_cscount)  # 处理其他情况（可选）
  ))


library(ggforce)
piedf <- as.data.frame.array(table(all_eGene_cscount_df_pie$all_eGene_cscount))
piedf$category <- row.names(piedf) 
colnames(piedf) <- c('value','category')

custom_colors <- c(
  "1" = "#0e3475",  # 深蓝色
  "2-10" = "#36528c",  # 深绿色
  "over_10" = "#99c8de"  # 深紫色
)

pie <- ggplot(piedf, aes(fill = category, amount = value)) +
  geom_arc_bar(stat = "pie", aes(x0 = 0, y0 = 0, r0 = 1, r = 2), color = "white") +
  coord_fixed()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(), 
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.title = element_blank(), 
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = 'none'
  ) +
  xlab("") + 
  ylab("") +
  scale_fill_manual(values = custom_colors) 

pdf(file = '/CIMA/Result/plot/trans_QTL/20250323_pie_transQTL.pdf',width = 6,height = 5)
print(pie)
dev.off()


library(ggplot2)

# 绘制条形图
bar <- ggplot(piedf, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity", color = "white", width = 0.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    legend.position = 'right'
  ) +
  xlab("") +
  ylab("Value") +
  scale_fill_manual(values = custom_colors)

# 保存为 PDF
pdf(file = '/CIMA/Result/plot/trans_QTL/20250323_bar_percentage_transQTL.pdf', width = 2, height = 8)
print(bar)
dev.off()

#画堆叠柱状图
CT_full_list <- read.csv('/CIMA/Result/summary/20250323_trans_eQTL_number_by_celltype.csv',row.names = 1)
plotting_sort <- read.csv('/CIMA/Data/69_celltype_sort.csv',row.names = 1)
L1_L4_mapping <- read.csv('/CIMA/Data/L1_L4.csv')
plotting_sort <- plotting_sort[plotting_sort$celltype %in% CT_full_list$V1,]


CT_full_list <-merge(CT_full_list,L1_L4_mapping,by.x = 'V1',by.y = 'celltype' )
CT_full_list$V1 <- factor(CT_full_list$V1,levels = plotting_sort)
CT_full_list$category =  factor(CT_full_list$category , levels = c("CD4 T", "CD8 T & unconvensional T", "NK&ILC", "B", "Myeloid",'HSPC'))
CT_full_list$diff = CT_full_list$eGene_num - CT_full_list$cs_eGene_num


barplot <- ggplot(CT_full_list, aes(x = V1, y = eGene_num)) +
  geom_bar(fill = "grey",stat = "identity", position = "stack") +
  geom_bar(aes(y = diff, fill = category), stat = "identity", position = "stack") +
  scale_fill_manual(values = c('#237b9f','#8ed1e5','#cca9cf','#f8b945','#91bd4a','#3ab37b')) +  # 使用tab20调色板
  #scale_x_discrete(labels = x_labels) +
  theme_classic() +
  labs(x = "", y = "Number of trans-eGenes", fill = "L1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())+
  scale_y_continuous(expand = c(0, 0.1))

pdf(file = '/CIMA/Result/plot/trans_QTL/20250323_barplot_transQTL.pdf',width = 10,height = 3)
print(barplot)
dev.off()

#画corr
eQTL_num = read.csv('/CIMA/Result/summary/20250212_eGene_number_csnumber_CT.csv',row.names = 1)
colnames(eQTL_num) <- paste0('cis_',colnames(eQTL_num))

CT_full_list <-merge(CT_full_list,eQTL_num,by.x = 'V1',by.y = 'cis_V1' )
lm_model <- lm(cis_eGene_num ~ eGene_num, data = CT_full_list )
# 提取相关性系数和p值
cor_coef <- round(cor(CT_full_list$cis_eGene_num, CT_full_list$eGene_num), 2)
p_value <- round(summary(lm_model)$coefficients[2, 4],30)

p <- ggplot(CT_full_list , aes(x = cis_eGene_num, y = eGene_num)) +
  geom_smooth(method = "lm",se = TRUE,show.legend = FALSE,color = "#466365",fill = "gray",linetype = "dashed",linewidth= 2)+
  labs(x = "cis_eGene_num", y = "trans_eGene_num", color = "L1") +
  theme_bw()+
  theme(panel.grid = element_blank())

# 添加回归直线和 95% 区间
p <- p + geom_point(aes(color = category),size = 3,show.legend = FALSE) + scale_color_manual(values = c("#237b9f"  ,"#8ed1e5","#cca9cf" ,"#f8b945" ,"#91bd4a","#3ab37b"))+theme(text = element_text(size = 12))+
  geom_text(x = min(CT_full_list$cis_eGene_num), y = max(CT_full_list$eGene_num)-3, 
            label = paste("Correlation:", cor_coef, "\nP-value:", p_value), 
            hjust = 0, vjust = 0, size = 5)

p
ggsave("/CIMA/Result/plot/trans_QTL/20250323_trans_cis_eQTL_corr.pdf", plot = p, width = 5, height = 5, dpi = 300, device = "pdf")

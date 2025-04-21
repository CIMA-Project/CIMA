library(ggplot2)
library(RColorBrewer)

eQTL_num = read.csv('/CIMA/Result/summary/20250212_eGene_number_csnumber_CT.csv',row.names = 1)
cell_percentage <- read.csv("/CIMA/Data/cell_proportions/proportions_final_annotation_RNA.csv")
mean_percentage <- cell_percentage %>%
  group_by(final_annotation) %>%
  summarize(mean_percentage = mean(proportion, na.rm = TRUE))

eQTL_num$cell_proportion <- 0
for (celltype in eQTL_num$V1){
  eQTL_num[eQTL_num$V1 == celltype,"cell_proportion"] <-  mean_percentage$mean_percentage[mean_percentage$final_annotation == celltype]
}

L1_L4 = read.csv('/CIMA/Data/L1_L4.csv')
eQTL_num$category <- 'B'
for (celltype in L1_L4$celltype){
  eQTL_num[eQTL_num$V1 == celltype,"category"] <-  L1_L4$category[L1_L4$celltype == celltype]
}

eQTL_num$category <- factor(eQTL_num$category , levels = c("CD4 T", "CD8 T & unconvensional T", "NK&ILC", "B", "Myeloid",'HSPC'))
# 计算线性回归模型
lm_model <- lm(eGene_num ~ cell_proportion, data = eQTL_num )
# 提取相关性系数和p值
cor_coef <- round(cor(eQTL_num$cell_proportion, eQTL_num$eGene_num), 2)
p_value <- round(summary(lm_model)$coefficients[2, 4],30)

p <- ggplot(eQTL_num , aes(x = cell_proportion, y = eGene_num)) +
  geom_smooth(method = "lm",se = TRUE,show.legend = FALSE,color = "#466365",fill = "gray",linetype = "dashed",linewidth= 2)+
  labs(x = "cell_proportion_RNA", y = "eGene_num", color = "L1") +
  theme_bw()+
  theme(panel.grid = element_blank())

# 添加回归直线和 95% 区间
p <- p + geom_point(aes(color = category),size = 3,show.legend = FALSE) + scale_color_manual(values = c("#237b9f"  ,"#8ed1e5","#cca9cf" ,"#f8b945" ,"#91bd4a","#3ab37b"))+theme(text = element_text(size = 12))+
  geom_text(x = min(eQTL_num$cell_proportion), y = max(eQTL_num$eGene_num)+1000, 
            label = paste("Correlation:", cor_coef, "\nP-value:", p_value), 
            hjust = 0, vjust = 0, size = 5)

p
ggsave("/CIMA/Result/plot/20250212_eQTL_cell_pro_corr.pdf", plot = p, width = 5, height = 5, dpi = 300, device = "pdf")


##处理caQTL
caQTL_num = read.csv('/CIMA/Result/summary/20250212_caPeak_number_csnumber_CT.csv',row.names = 1)
cell_percentage <- read.csv("/CIMA/Data/cell_proportions/proportions_final_annotation_ATAC.csv")
mean_percentage <- cell_percentage %>%
  group_by(final_annotation) %>%
  summarize(mean_percentage = mean(proportion, na.rm = TRUE))

caQTL_num$cell_proportion <- 0
for (celltype in caQTL_num$V1){
  caQTL_num[caQTL_num$V1 == celltype,"cell_proportion"] <-  mean_percentage$mean_percentage[mean_percentage$final_annotation == celltype]
}

L1_L4 = read.csv('/CIMA/Data/L1_L4.csv')
caQTL_num$category <- 'B'
for (celltype in L1_L4$celltype){
  caQTL_num[caQTL_num$V1 == celltype,"category"] <-  L1_L4$category[L1_L4$celltype == celltype]
}

caQTL_num$category <- factor(caQTL_num$category , levels = c("CD4 T", "CD8 T & unconvensional T", "NK&ILC", "B", "Myeloid"))
# 计算线性回归模型
lm_model <- lm(caPeak_num ~ cell_proportion, data = caQTL_num )
# 提取相关性系数和p值
cor_coef <- round(cor(caQTL_num$cell_proportion, caQTL_num$caPeak_num), 2)
p_value <- round(summary(lm_model)$coefficients[2, 4],30)

p <- ggplot(caQTL_num , aes(x = cell_proportion, y = caPeak_num)) +
  geom_smooth(method = "lm",se = TRUE,show.legend = FALSE,color = "#466365",fill = "gray",linetype = "dashed",linewidth= 2)+
  labs(x = "cell_proportion_ATAC", y = "caPeak_num", color = "L1") +
  theme_bw()+
  theme(panel.grid = element_blank())

# 添加回归直线和 95% 区间
p <- p + geom_point(aes(color = category),size = 3,show.legend = FALSE)  + scale_color_manual(values = c("#237b9f"  ,"#8ed1e5","#cca9cf" ,"#f8b945" ,"#91bd4a"))+theme(text = element_text(size = 12))+
  geom_text(x = min(caQTL_num$cell_proportion), y = max(caQTL_num$caPeak_num-5000), 
            label = paste("Correlation:", cor_coef, "\nP-value:", p_value), 
            hjust = 0, vjust = 0, size = 5)

p
ggsave("/CIMA/Result/plot/20250212_caQTL_cell_pro_corr.pdf", plot = p, width = 5, height = 5, dpi = 300, device = "pdf")


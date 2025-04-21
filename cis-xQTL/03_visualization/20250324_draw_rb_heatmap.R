library(pheatmap)
library(viridis)
library(ggplot2)

custom_colors <- colorRampPalette(c("#FFC787", "white", "#540D6E"))(100)

eGene <- read.csv("/CIMA/Result/summary/20250212_eGene_number_csnumber_CT.csv",row.names = 1)
eGene  <- eGene [eGene$eGene_num > 100,]

celltype_sort = read.csv('/CIMA/Data/69_celltype_sort.csv')
celltype_sort = celltype_sort[celltype_sort$celltype %in% eGene$V1,]
eQTL_rb_csv = read.csv('/CIMA/Result/rb_test/20250324_eQTL_rb.csv',row.names = 1)

row.names(eQTL_rb_csv) <- gsub('-','.',row.names(eQTL_rb_csv))
for (CT in row.names(eQTL_rb_csv)){
  eQTL_rb_csv[CT,CT] <- 1
}

eQTL_rb_csv <- eQTL_rb_csv[gsub("-", ".", celltype_sort$celltype),gsub("-", ".", celltype_sort$celltype)] 
pheatmap(eQTL_rb_csv,cluster_rows = F,cluster_cols = F,color = custom_colors)

min(eQTL_rb_csv)
# 提取数据范围
min_val <- min(eQTL_rb_csv)
max_val <- max(eQTL_rb_csv)

# 生成breaks
breaks_1 <- seq(min_val, 0.8, length.out = 51)  # 50个区间，51个break点
breaks_2 <- seq(0.8, max_val, length.out = 51)

# 合并breaks，去重以避免重复0.8
breaks <- c(breaks_1, breaks_2[-1])

pheatmap(eQTL_rb_csv,cluster_rows = F,cluster_cols = F,color = custom_colors,breaks = breaks)
a1 <- pheatmap(eQTL_rb_csv,cluster_rows = F,cluster_cols = F,color = custom_colors,breaks = breaks)
ggsave(a1,filename = '/CIMA/Result/plot/rb/20250304_eQTL_rb_61ct.pdf',width = 10,height = 10)


#画caQTL
caPeak <- read.csv("/CIMA/Result/summary/20250212_caPeak_number_csnumber_CT.csv",row.names = 1)
caPeak <- caPeak[caPeak$caPeak_num > 100,]

celltype_sort = read.csv('/CIMA/Data/69_celltype_sort.csv')
celltype_sort = celltype_sort[celltype_sort$celltype %in% caPeak$V1,]
eQTL_rb_csv = read.csv('/CIMA/Result/rb_test/20250324_caQTL_rb.csv',row.names = 1)


row.names(eQTL_rb_csv) <- gsub('-','.',row.names(eQTL_rb_csv))
for (CT in row.names(eQTL_rb_csv)){
  eQTL_rb_csv[CT,CT] <- 1
}


eQTL_rb_csv <- eQTL_rb_csv[gsub("-", ".", celltype_sort$celltype),gsub("-", ".", celltype_sort$celltype)] 
pheatmap(eQTL_rb_csv,cluster_rows = F,cluster_cols = F,color = custom_colors)

min(eQTL_rb_csv)
# 提取数据范围
min_val <- min(eQTL_rb_csv)
max_val <- max(eQTL_rb_csv)

# 生成breaks
breaks_1 <- seq(min_val, 0.8, length.out = 51)  # 50个区间，51个break点
breaks_2 <- seq(0.8, max_val, length.out = 51)

# 合并breaks，去重以避免重复0.8
breaks <- c(breaks_1, breaks_2[-1])

pheatmap(eQTL_rb_csv,cluster_rows = F,cluster_cols = F,color = custom_colors,breaks = breaks)
a1 <- pheatmap(eQTL_rb_csv,cluster_rows = F,cluster_cols = F,color = custom_colors,breaks = breaks)
ggsave(a1,filename = '/CIMA/Result/plot/rb/20250304_caQTL_rb_35ct.pdf',width = 10,height = 10)

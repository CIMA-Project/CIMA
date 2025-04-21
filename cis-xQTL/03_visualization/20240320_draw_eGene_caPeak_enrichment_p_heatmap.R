library(pheatmap)
library(RColorBrewer)

celltype_sort = read.csv('/CIMA/Data/69_celltype_sort.csv')
eQTL_enrich_p_df <- read.csv('/CIMA/Result/summary/20250206_eGene_enrich_p_df.csv',row.names = 1)
#eQTL_enrich_p_df <- apply(eQTL_enrich_p_df, c(1, 2), function(x) ifelse(!is.na(x), -log10(x), NA))
#eQTL_enrich_p_df <- as.data.frame(t(eQTL_enrich_p_df))
eQTL_enrich_p_df <- eQTL_enrich_p_df[celltype_sort$celltype,gsub("-", ".", celltype_sort$celltype)] 

sum(eQTL_enrich_p_df < 0.05, na.rm = TRUE)

# 使用 RdYlBu 调色板
#colors <- colorRampPalette(brewer.pal(9, "BuPu"))(100)
colors <- colorRampPalette(c("#F7FCFD","#F7FCFD","#F7FCFD","#F7FCFD","#F7FCFD","#9ECAE1","#6BAED6","#4292C6","#2171B5","#08519C"))(100)
# 设置断点，1.3是转折点，确保没有重复值
min_val <- min(eQTL_enrich_p_df, na.rm = TRUE)
max_val <- max(eQTL_enrich_p_df, na.rm = TRUE)

# 创建分段：从 min 到 1.3，再从 1.3 到 max
breaks <- c(seq(min_val, 0.05, length.out = 50), 
            seq(0.05 + .Machine$double.eps, max_val, length.out = 51))

pdf('/CIMA/Result/plot/basic_number/20250320_eGene_enrichment_P.pdf',width = 12,height = 10)
# 创建一个pheatmap
pheatmap(eQTL_enrich_p_df, 
         color = colors, 
         breaks = breaks, 
         na_col = "grey",  # NA值的颜色是白色，或者你可以选择“transparent”
         legend = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         legend_breaks = c(min_val, 0.05,0.5, max_val),
         border_color = NA)
dev.off()

###

celltype_sort = read.csv('/CIMA/Data/69_celltype_sort.csv')
caQTL_enrich_p_df <- read.csv('/CIMA/Result/summary/20250206_caPeak_enrich_p_df.csv',row.names = 1)
#caQTL_enrich_p_df <- apply(caQTL_enrich_p_df, c(1, 2), function(x) ifelse(!is.na(x), -log10(x), NA))
#caQTL_enrich_p_df <- as.data.frame(t(caQTL_enrich_p_df))
celltype_sort <- celltype_sort[celltype_sort$celltype %in% row.names(caQTL_enrich_p_df),]
caQTL_enrich_p_df <- caQTL_enrich_p_df[celltype_sort$celltype,gsub("-", ".", celltype_sort$celltype)] 

sum(caQTL_enrich_p_df < 0.05, na.rm = TRUE)

# 使用 RdYlBu 调色板
#colors <- colorRampPalette(brewer.pal(9, "BuPu"))(100)
colors <- colorRampPalette(c("#F7FCFD","#F7FCFD","#F7FCFD","#F7FCFD","#F7FCFD","#9ECAE1","#6BAED6","#4292C6","#2171B5","#08519C"))(100)
# 设置断点，1.3是转折点，确保没有重复值
min_val <- min(caQTL_enrich_p_df, na.rm = TRUE)
max_val <- max(caQTL_enrich_p_df, na.rm = TRUE)

# 创建分段：从 min 到 1.3，再从 1.3 到 max
breaks <- c(seq(min_val, 0.05, length.out = 50), 
            seq(0.05 + .Machine$double.eps, max_val, length.out = 51))

pdf('/CIMA/Result/plot/basic_number/20250320_caPeak_enrichment_P.pdf',width = 12,height = 10)
# 创建一个pheatmap
pheatmap(caQTL_enrich_p_df, 
         color = colors, 
         breaks = breaks, 
         na_col = "grey",  # NA值的颜色是白色，或者你可以选择“transparent”
         legend = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         legend_breaks = c(min_val, 0.05,0.5, max_val),
         border_color = NA)
dev.off()


library(pheatmap)
library(RColorBrewer)

plot_df <- read.csv('/CIMA/Result/summary/20250320_celltype_independent_percentage_eQTL.csv',row.names = 1)
min(plot_df,na.rm = T)
max(plot_df,na.rm = T)
mean(as.matrix(plot_df),na.rm = T)

# Define the breaks
# 50 intervals between 0 and 0.7
colors = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
low_breaks <- seq(0, 0.7, length.out = 51)

# 50 intervals between 0.7 and 1
high_breaks <- seq(0.7, 1, length.out = 51)

# Combine both breaks
breaks <- c(low_breaks, high_breaks[-1])  # Remove the duplicated 0.7 at the transition

pdf('/CIMA/Result/plot/basic_number/20250320_celltype_independent_percentage_eQTL.pdf',width = 10,height = 10)
# 创建一个pheatmap
pheatmap(plot_df, 
         color = colors, 
         na_col = "grey",  # NA值的颜色是白色，或者你可以选择“transparent”
         legend = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         breaks = breaks,
         border_color = NA)
dev.off()

plot_df <- read.csv('/CIMA/Result/summary/20250320_celltype_independent_percentage_caQTL.csv',row.names = 1)
min(plot_df,na.rm = T)
max(plot_df,na.rm = T)
mean(as.matrix(plot_df),na.rm = T)
# Define the breaks
# 50 intervals between 0 and 0.7
colors = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
low_breaks <- seq(0, 0.7, length.out = 51)

# 50 intervals between 0.7 and 1
high_breaks <- seq(0.7, 1, length.out = 51)

# Combine both breaks
breaks <- c(low_breaks, high_breaks[-1])  # Remove the duplicated 0.7 at the transition

pdf('/CIMA/Result/plot/basic_number/20250320_celltype_independent_percentage_caQTL.pdf',width = 10,height = 10)
# 创建一个pheatmap
pheatmap(plot_df, 
         color = colors, 
         na_col = "grey",  # NA值的颜色是白色，或者你可以选择“transparent”
         legend = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         breaks = breaks,
         border_color = NA)
dev.off()


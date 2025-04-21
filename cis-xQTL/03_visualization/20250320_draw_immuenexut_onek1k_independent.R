library(reshape2)
library(pheatmap)

plot_df <- read.csv('/CIMA/Result/summary/20250320_onek1k_immunexut_independent_percentage_eQTL.csv',row.names = 1)
plot_df <- as.data.frame(t(plot_df))

plot_df <- plot_df[,c('CD4_T','CD8_T',"NK","B","Myeloid")]
breaks <- seq(0.6, 0.9, length.out = 100)
pdf('/CIMA/Result/plot/20250320_immuenexut_onek1k_vs_CIMA_independent_percentage.pdf',width = 6,height = 3)
# 绘制热图
pheatmap(
  plot_df,
  display_numbers = TRUE,  # 显示值
  cluster_rows = FALSE,               # 禁用行聚类
  cluster_cols = FALSE,               # 禁用列聚类
  breaks = breaks,                    # 设置值域
  #color = colorRampPalette(c("white", "red"))(100),  # 设置颜色
  fontsize_number = 10 )               # 设置值的大小
dev.off()

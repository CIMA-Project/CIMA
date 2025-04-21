library(reshape2)
library(pheatmap)

pi1_1 <- read.csv('/CIMA/Result/pi1_test/public/20250219_onek1k_pi0.csv')
pi1_2 <- read.csv('/CIMA/Result/pi1_test/public/20250219_immuenexut_pi0.csv')

pi1_1$V1 <- 1 - pi1_1$V1
pi1_2$V1 <- 1 - pi1_2$V1

pi1_1$celltype.s <- factor(pi1_1$celltype.s,levels = c('CD4_T','CD8_T','NK','B','Myeloid'))
pi1_1$celltype.compared <- factor(pi1_1$celltype.compare,levels = c('CD4_T_own','CD8_T_own','NK_own','B_own','Myeloid_own'))

df_wide <- dcast(pi1_1, celltype.s ~ celltype.compared, value.var = "V1")
rownames(df_wide) <- df_wide$celltype.s
df_wide$celltype.s <- NULL

# 定义值域（从 0.5 到 1）
breaks <- seq(0.5, 0.85, length.out = 100)
pdf('/CIMA/Result/plot/20250320_onek1k_vs_CIMA_pi1.pdf',width = 6,height = 5)
# 绘制热图
pheatmap(
  df_wide,
  display_numbers = TRUE,  # 显示值
  cluster_rows = FALSE,               # 禁用行聚类
  cluster_cols = FALSE,               # 禁用列聚类
  breaks = breaks,                    # 设置值域
  #color = colorRampPalette(c("white", "red"))(100),  # 设置颜色
  fontsize_number = 10 )               # 设置值的大小
dev.off()


pi1_2$celltype.s <- factor(pi1_2$celltype.s,levels = c('CD4_T','CD8_T','NK','B','Myeloid'))
pi1_2$celltype.compared <- factor(pi1_2$celltype.compare,levels = c('CD4_T_own','CD8_T_own','NK_own','B_own','Myeloid_own'))

df_wide <- dcast(pi1_2, celltype.s ~ celltype.compared, value.var = "V1")
rownames(df_wide) <- df_wide$celltype.s
df_wide$celltype.s <- NULL

# 定义值域（从 0.5 到 1）
breaks <- seq(0.5, 0.85, length.out = 100)
pdf('/CIMA/Result/plot/20250320_immuenexut_vs_CIMA_pi1.pdf',width = 6,height = 5)
# 绘制热图
pheatmap(
  df_wide,
  display_numbers = TRUE,  # 显示值
  cluster_rows = FALSE,               # 禁用行聚类
  cluster_cols = FALSE,               # 禁用列聚类
  breaks = breaks,                    # 设置值域
  #color = colorRampPalette(c("white", "red"))(100),  # 设置颜色
  fontsize_number = 10 )               # 设置值的大小
dev.off()

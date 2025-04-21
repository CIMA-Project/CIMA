library(ggplot2)
library(gridExtra)
library(dplyr)

# 为 celltype_pair 定义颜色映射
color_mapping <- c(
  "CD4_T_xxx_CD4_T_own" = "#237b9f",
  "CD8_T_xxx_CD8_T_own" = "#8ed1e5",
  "NK_xxx_NK_own" = "#cca9cf",
  "Myeloid_xxx_Myeloid_own" = "#91bd4a",
  "B_xxx_B_own" = "#f8b945"
)

onek1k_vs_own = read.csv('/CIMA/Result/eQTL_L1_downstream/20250219_onek1k_vs_own.csv')
onek1k_vs_own = onek1k_vs_own[onek1k_vs_own$celltype_pair %in% c('CD4_T_xxx_CD4_T_own','CD8_T_xxx_CD8_T_own','NK_xxx_NK_own','Myeloid_xxx_Myeloid_own','B_xxx_B_own'),]
onek1k_vs_own = onek1k_vs_own[,c('celltype_pair','slope','beta_public')]
colnames(onek1k_vs_own) = c('celltype_pair','beta_CIMA','beta_onek1k')

# 计算相关性系数，并绘制图形
plots <- lapply(unique(onek1k_vs_own$celltype_pair), function(cat) {
  sub_df <- subset(onek1k_vs_own, celltype_pair == cat)
  
  # 计算相关性系数
  cor_val <- round(cor(sub_df$beta_CIMA, sub_df$beta_onek1k), 2)
  # 创建散点图
  p <- ggplot(sub_df, aes(x = beta_CIMA, y = beta_onek1k)) +
    geom_point(aes(color = celltype_pair), size = 1) +  # 按 celltype_pair 上色
    ggtitle(paste("Category", cat, "\nCorrelation:", cor_val)) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),  # 去除网格线
      plot.title = element_text(hjust = 0.5),  # 标题居中
      legend.position = "none"  # 去除图例
    ) +
    scale_color_manual(values = color_mapping) +  # 按照 celltype_pair 设置颜色
    geom_vline(xintercept = 0, color = "black", linetype = "dashed") +  # 添加 vertical line 在 beta_CIMA = 0
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") + # 添加 horizontal line 在 beta_onek1k = 0
    xlim(-2, 2)+
    ylim(-1, 1)
  return(p)
})

# 将五张图拼成一行
p_hold <- grid.arrange(grobs = plots, ncol = 5)
dev.off()

pdf("/CIMA/Result/plot/20250220_onek1k_vs_cima_beta.pdf", width = 15, height = 4)  # 你可以调整宽度和高度
grid.arrange(grobs = plots, ncol = 5)
# 关闭设备，保存文件
dev.off()


ImmuNexUT_vs_own = read.csv('/CIMA/Result/eQTL_L1_downstream/20250219_immuenexut_vs_own.csv')
ImmuNexUT_vs_own = ImmuNexUT_vs_own[ImmuNexUT_vs_own$celltype_pair %in% c('CD4_T_xxx_CD4_T_own','CD8_T_xxx_CD8_T_own','NK_xxx_NK_own','Myeloid_xxx_Myeloid_own','B_xxx_B_own'),]
ImmuNexUT_vs_own = ImmuNexUT_vs_own[,c('celltype_pair','slope','beta_public')]
colnames(ImmuNexUT_vs_own) = c('celltype_pair','beta_CIMA','beta_ImmuNexUT')
# 计算相关性系数，并绘制图形
# 计算相关性系数，并绘制图形
plots <- lapply(unique(ImmuNexUT_vs_own$celltype_pair), function(cat) {
  sub_df <- subset(ImmuNexUT_vs_own, celltype_pair == cat)
  
  # 使用 cor.test 计算斯皮尔曼相关系数
  cor_val <- round(cor(sub_df$beta_CIMA, sub_df$beta_ImmuNexUT), 2)
  
  # 创建散点图
  p <- ggplot(sub_df, aes(x = beta_CIMA, y = beta_ImmuNexUT)) +
    geom_point(aes(color = celltype_pair), size = 1) +  # 按 celltype_pair 上色
    ggtitle(paste("Category", cat, "\nCorrelation:", cor_val)) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),  # 去除网格线
      plot.title = element_text(hjust = 0.5),  # 标题居中
      legend.position = "none"  # 去除图例
    ) +
    scale_color_manual(values = color_mapping) +  # 按照 celltype_pair 设置颜色
    geom_vline(xintercept = 0, color = "black", linetype = "dashed") +  # 添加 vertical line 在 beta_CIMA = 0
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") + # 添加 horizontal line 在 beta_onek1k = 0
    xlim(-2, 2)+
    ylim(-2, 2)
  return(p)
})


# 将五张图拼成一行
p_hold <- grid.arrange(grobs = plots, ncol = 5)
dev.off()

pdf("/CIMA/Result/plot/20250220_immuenexut_vs_cima_beta.pdf", width = 15, height = 4)  # 你可以调整宽度和高度
grid.arrange(grobs = plots, ncol = 5)
# 关闭设备，保存文件
dev.off()


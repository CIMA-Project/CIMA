setwd('/CIMA/Script/')
library(tidyverse)
library(dplyr)
library(bestNormalize)
library(parallel)

# 假设你已经加载了所需的数据和 orderNorm_gene 函数
orderNorm_gene <- function(gene.name, expr.mat)
{
  expr_norm <- expr.mat[,gene.name] %>% orderNorm()
  df <- data.frame(gene = expr_norm$x.t)
  colnames(df) <- gene.name  
  return(df)
}

celltype_eQTL = c("B","Monocyte")

for (celltype in celltype_eQTL) {
  
  # 根据细胞类型构建第一个文件路径（假设文件路径为 '/pseudobulk' 文件夹）
  pseudobulk_matrix_path1 <- paste0('/CIMA/Data/dynamic/pseudobulk/pseudobulk/', celltype, '.csv')
  
  # 读取第一个文件
  pseudobulk_matrix_1 <- read.csv(pseudobulk_matrix_path1, row.names = 1)
  
  # 使用 mclapply() 进行并行化处理
  norm_list1 <- mclapply(colnames(pseudobulk_matrix_1), orderNorm_gene, expr.mat = pseudobulk_matrix_1, mc.cores = detectCores())
  
  # 将处理后的结果绑定成一个数据框
  expr1 <- bind_cols(norm_list1)
  row.names(expr1) <- row.names(pseudobulk_matrix_1)
  
  # 输出到第一个文件夹（假设路径为 '/normal_dis'）
  output_path1 <- paste0('/CIMA/Data/dynamic/pseudobulk/normal_dis/', celltype, '.csv')
  write.csv(expr1, file = output_path1)
  
  # 根据细胞类型构建第二个文件路径（假设是 '/top2000_pseudobulk' 文件夹）
  pseudobulk_matrix_path2 <- paste0('/CIMA/Data/dynamic/pseudobulk/top2000_pseudobulk/', celltype, '.csv')
  
  # 读取第二个文件
  pseudobulk_matrix_2 <- read.csv(pseudobulk_matrix_path2, row.names = 1)
  
  # 使用 mclapply() 进行并行化处理
  norm_list2 <- mclapply(colnames(pseudobulk_matrix_2), orderNorm_gene, expr.mat = pseudobulk_matrix_2, mc.cores = detectCores())
  
  # 将处理后的结果绑定成一个数据框
  expr2 <- bind_cols(norm_list2)
  row.names(expr2) <- row.names(pseudobulk_matrix_2)
  
  # 输出到第二个文件夹（假设路径为 '/normal_dis_top2000'）
  output_path2 <- paste0('/CIMA/Data/dynamic/pseudobulk/top2000_normal_dis/', celltype, '.csv')
  write.csv(expr2, file = output_path2)
  
  # 打印日志
  cat("Processed and saved for celltype:", celltype, "\n")
}

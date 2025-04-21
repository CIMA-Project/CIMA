setwd('/CIMA/Script/')
library(tidyverse)
library(dplyr)
library(bestNormalize)
library(parallel)

# 鍋囪浣犲凡缁忓姞杞戒簡鎵€闇€鐨勬暟鎹拰 orderNorm_gene 鍑芥�?
orderNorm_gene <- function(gene.name, expr.mat)
{
  expr_norm <- expr.mat[,gene.name] %>% orderNorm()
  df <- data.frame(gene = expr_norm$x.t)
  colnames(df) <- gene.name  
  return(df)
}

celltype_eQTL = read.csv('/CIMA/Data/20241230_xQTL_cell_sample_number/20241230_celltype_70_scRNA.csv',row.names = )

for (celltype in celltype_eQTL$final_annotation) {
  
  # 鏍规嵁缁嗚優绫诲瀷鏋勫缓绗竴涓枃浠惰矾寰勶紙鍋囪鏂囦欢璺緞涓? '/pseudobulk' 鏂囦欢澶癸級
  pseudobulk_matrix_path1 <- paste0('/CIMA/Data/eQTL/pseudobulk/', celltype, '.csv')
  
  # 璇诲彇绗竴涓枃浠?
  pseudobulk_matrix_1 <- read.csv(pseudobulk_matrix_path1, row.names = 1)
  
  # 浣跨�? mclapply() 杩涜骞惰鍖栧鐞?
  norm_list1 <- mclapply(colnames(pseudobulk_matrix_1), orderNorm_gene, expr.mat = pseudobulk_matrix_1, mc.cores = detectCores())
  
  # 灏嗗鐞嗗悗鐨勭粨鏋滅粦瀹氭垚涓€涓暟鎹�?
  expr1 <- bind_cols(norm_list1)
  row.names(expr1) <- row.names(pseudobulk_matrix_1)
  
  # 杈撳嚭鍒扮涓€涓枃浠跺す锛堝亣璁捐矾寰勪�? '/normal_dis'�??
  output_path1 <- paste0('/CIMA/Data/eQTL/normal_dis/', celltype, '.csv')
  write.csv(expr1, file = output_path1)
  
  # 鏍规嵁缁嗚優绫诲瀷鏋勫缓绗簩涓枃浠惰矾寰勶紙鍋囪鏄? '/top2000_pseudobulk' 鏂囦欢澶癸級
  pseudobulk_matrix_path2 <- paste0('/CIMA/Data/eQTL/top2000_pseudobulk/', celltype, '.csv')
  
  # 璇诲彇绗簩涓枃浠?
  pseudobulk_matrix_2 <- read.csv(pseudobulk_matrix_path2, row.names = 1)
  
  # 浣跨�? mclapply() 杩涜骞惰鍖栧鐞?
  norm_list2 <- mclapply(colnames(pseudobulk_matrix_2), orderNorm_gene, expr.mat = pseudobulk_matrix_2, mc.cores = detectCores())
  
  # 灏嗗鐞嗗悗鐨勭粨鏋滅粦瀹氭垚涓€涓暟鎹�?
  expr2 <- bind_cols(norm_list2)
  row.names(expr2) <- row.names(pseudobulk_matrix_2)
  
  # 杈撳嚭鍒扮浜屼釜鏂囦欢澶癸紙鍋囪璺緞涓? '/normal_dis_top2000'�??
  output_path2 <- paste0('/CIMA/Data/eQTL/top2000_normal_dis/', celltype, '.csv')
  write.csv(expr2, file = output_path2)
  
  # 鎵撳嵃鏃ュ織
  cat("Processed and saved for celltype:", celltype, "\n")
}

celltype_caQTL = read.csv('/CIMA/Data/20241230_xQTL_cell_sample_number/20241230_celltype_70_scATAC.csv',row.names = )

for (celltype in celltype_caQTL$final_annotation) {
  
  # 鏍规嵁缁嗚優绫诲瀷鏋勫缓绗竴涓枃浠惰矾寰勶紙鍋囪鏂囦欢璺緞涓? '/pseudobulk' 鏂囦欢澶癸級
  pseudobulk_matrix_path1 <- paste0('/CIMA/Data/caQTL/pseudobulk/', celltype, '.csv')
  
  # 璇诲彇绗竴涓枃浠?
  pseudobulk_matrix_1 <- read.csv(pseudobulk_matrix_path1, row.names = 1)
  
  # 浣跨�? mclapply() 杩涜骞惰鍖栧鐞?
  norm_list1 <- mclapply(colnames(pseudobulk_matrix_1), orderNorm_gene, expr.mat = pseudobulk_matrix_1, mc.cores = detectCores())
  
  # 灏嗗鐞嗗悗鐨勭粨鏋滅粦瀹氭垚涓€涓暟鎹�?
  expr1 <- bind_cols(norm_list1)
  row.names(expr1) <- row.names(pseudobulk_matrix_1)
  
  # 杈撳嚭鍒扮涓€涓枃浠跺す锛堝亣璁捐矾寰勪�? '/normal_dis'�??
  output_path1 <- paste0('/CIMA/Data/caQTL/normal_dis/', celltype, '.csv')
  write.csv(expr1, file = output_path1)
  
  # 鏍规嵁缁嗚優绫诲瀷鏋勫缓绗簩涓枃浠惰矾寰勶紙鍋囪鏄? '/top2000_pseudobulk' 鏂囦欢澶癸級
  pseudobulk_matrix_path2 <- paste0('/CIMA/Data/caQTL/top2000_pseudobulk/', celltype, '.csv')
  
  # 璇诲彇绗簩涓枃浠?
  pseudobulk_matrix_2 <- read.csv(pseudobulk_matrix_path2, row.names = 1)
  
  # 浣跨�? mclapply() 杩涜骞惰鍖栧鐞?
  norm_list2 <- mclapply(colnames(pseudobulk_matrix_2), orderNorm_gene, expr.mat = pseudobulk_matrix_2, mc.cores = detectCores())
  
  # 灏嗗鐞嗗悗鐨勭粨鏋滅粦瀹氭垚涓€涓暟鎹�?
  expr2 <- bind_cols(norm_list2)
  row.names(expr2) <- row.names(pseudobulk_matrix_2)
  
  # 杈撳嚭鍒扮浜屼釜鏂囦欢澶癸紙鍋囪璺緞涓? '/normal_dis_top2000'�??
  output_path2 <- paste0('/CIMA/Data/caQTL/top2000_normal_dis/', celltype, '.csv')
  write.csv(expr2, file = output_path2)
  
  # 鎵撳嵃鏃ュ織
  cat("Processed and saved for celltype:", celltype, "\n")
}


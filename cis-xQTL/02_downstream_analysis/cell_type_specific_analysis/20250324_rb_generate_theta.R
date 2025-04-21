library(data.table)
library(glue)
library(parallel)

#计算eQTL
CT_for_rb_eQTL = read.csv('/CIMA/Data/20241230_xQTL_cell_sample_number/20241230_celltype_70_scRNA.csv')
CT_for_rb_eQTL = CT_for_rb_eQTL$final_annotation

#all_combinations <- combn(CT_for_rb_eQTL , 2, simplify = FALSE)
sample_overlap_df <- data.frame(matrix(NA, nrow = length(CT_for_rb_eQTL), ncol = length(CT_for_rb_eQTL)))
colnames(sample_overlap_df) <- CT_for_rb_eQTL
rownames(sample_overlap_df) <- CT_for_rb_eQTL

sample_number_df <- sample_overlap_df[,'CD4_Tn_CCR7',drop = FALSE]
colnames(sample_number_df) <- c('sample_number')

corr_list <- list()

# 读取所有数据框到一个列表中
expr_data_list <- lapply(CT_for_rb_eQTL, function(CT_ref) {
  data.frame(fread(glue('/CIMA/Data/eQTL/normal_dis/{CT_ref}.csv')), row.names = 1)
})
names(expr_data_list) <- CT_for_rb_eQTL

# 更新 sample_number_df 中的样本数量
sample_number_df$sample_number <- sapply(expr_data_list, nrow)
min(sample_number_df$sample_number)

for (CT_ref in CT_for_rb_eQTL){
  expr_CT_ref <- expr_data_list[[CT_ref]]
  for (CT_query in CT_for_rb_eQTL[CT_for_rb_eQTL != CT_ref]){
    #print(glue('ref: {CT_ref};query: {CT_query}'))
    expr_CT_query <- expr_data_list[[CT_query]]
    sample_overlap <- length(intersect(row.names(expr_CT_ref),row.names(expr_CT_query)))/sqrt(nrow(expr_CT_ref)*nrow(expr_CT_query))
    sample_overlap_df[match(CT_ref,row.names(sample_overlap_df)),CT_query] <- sample_overlap
    
  }
}

calculate_correlations <- function(pair) {
  CT_ref <- pair[1]
  CT_query <- pair[2]
  
  expr_CT_ref <- expr_data_list[[CT_ref]]
  expr_CT_query <- expr_data_list[[CT_query]]
  
  # 裁剪两个数据框，使其只保留重叠的行和共同的列
  common_rows <- intersect(row.names(expr_CT_ref), row.names(expr_CT_query))
  common_cols <- intersect(colnames(expr_CT_ref), colnames(expr_CT_query))
  
  expr_CT_ref_trimmed <- expr_CT_ref[common_rows, common_cols, drop=FALSE]
  expr_CT_query_trimmed <- expr_CT_query[common_rows, common_cols, drop=FALSE]
  
  # 计算每一列的皮尔逊相关系�?
  correlations <- sapply(common_cols, function(colname) {
    cor(expr_CT_ref_trimmed[[colname]], expr_CT_query_trimmed[[colname]])
  })
  
  # 返回样本重叠和相关系�?
  return(correlations)
}

pairwise_combinations <- combn(CT_for_rb_eQTL, 2)

# 获取核心�?
num_cores <- round(detectCores()/3,0)

# 使用 mclapply 进行并行计算 (适用于Unix类操作系�?)
results <- mclapply(1:ncol(pairwise_combinations), function(i) calculate_correlations(pairwise_combinations[, i]), mc.cores = num_cores)

pair <- c()
for (i in 1:ncol(pairwise_combinations)){ 
  pair_tmp <- paste0(pairwise_combinations[1, i],'*',pairwise_combinations[2, i])
  pair <- c(pair,pair_tmp)
}

names(results) <- pair
saveRDS(results, "/CIMA/Result/rb_test/20250324_eQTL_rb_corr.rds")
write.csv(sample_number_df, "/CIMA/Result/rb_test/20250324_eQTL_sample_number.csv")
write.csv(sample_overlap_df, "/CIMA/Result/rb_test/20250324_eQTL_sample_overlap.csv")

#计算caQTL
CT_for_rb_eQTL = read.csv('/CIMA/Data/20241230_xQTL_cell_sample_number/20241230_celltype_70_scATAC.csv')
CT_for_rb_eQTL = CT_for_rb_eQTL$final_annotation

#all_combinations <- combn(CT_for_rb_eQTL , 2, simplify = FALSE)
sample_overlap_df <- data.frame(matrix(NA, nrow = length(CT_for_rb_eQTL), ncol = length(CT_for_rb_eQTL)))
colnames(sample_overlap_df) <- CT_for_rb_eQTL
rownames(sample_overlap_df) <- CT_for_rb_eQTL

sample_number_df <- sample_overlap_df[,'CD4_Tn_CCR7',drop = FALSE]
colnames(sample_number_df) <- c('sample_number')

corr_list <- list()

# 读取所有数据框到一个列表中
expr_data_list <- lapply(CT_for_rb_eQTL, function(CT_ref) {
  data.frame(fread(glue('/CIMA/Data/caQTL/normal_dis/{CT_ref}.csv')), row.names = 1)
})
names(expr_data_list) <- CT_for_rb_eQTL

# 更新 sample_number_df 中的样本数量
sample_number_df$sample_number <- sapply(expr_data_list, nrow)
min(sample_number_df$sample_number)

for (CT_ref in CT_for_rb_eQTL){
  expr_CT_ref <- expr_data_list[[CT_ref]]
  for (CT_query in CT_for_rb_eQTL[CT_for_rb_eQTL != CT_ref]){
    #print(glue('ref: {CT_ref};query: {CT_query}'))
    expr_CT_query <- expr_data_list[[CT_query]]
    sample_overlap <- length(intersect(row.names(expr_CT_ref),row.names(expr_CT_query)))/sqrt(nrow(expr_CT_ref)*nrow(expr_CT_query))
    sample_overlap_df[match(CT_ref,row.names(sample_overlap_df)),CT_query] <- sample_overlap
    
  }
}

calculate_correlations <- function(pair) {
  CT_ref <- pair[1]
  CT_query <- pair[2]
  
  expr_CT_ref <- expr_data_list[[CT_ref]]
  expr_CT_query <- expr_data_list[[CT_query]]
  
  # 裁剪两个数据框，使其只保留重叠的行和共同的列
  common_rows <- intersect(row.names(expr_CT_ref), row.names(expr_CT_query))
  common_cols <- intersect(colnames(expr_CT_ref), colnames(expr_CT_query))
  
  expr_CT_ref_trimmed <- expr_CT_ref[common_rows, common_cols, drop=FALSE]
  expr_CT_query_trimmed <- expr_CT_query[common_rows, common_cols, drop=FALSE]
  
  # 计算每一列的皮尔逊相关系�?
  correlations <- sapply(common_cols, function(colname) {
    cor(expr_CT_ref_trimmed[[colname]], expr_CT_query_trimmed[[colname]])
  })
  
  # 返回样本重叠和相关系�?
  return(correlations)
}

pairwise_combinations <- combn(CT_for_rb_eQTL, 2)

# 获取核心�?
num_cores <- round(detectCores()/3,0)

# 使用 mclapply 进行并行计算 (适用于Unix类操作系�?)
results <- mclapply(1:ncol(pairwise_combinations), function(i) calculate_correlations(pairwise_combinations[, i]), mc.cores = num_cores)

pair <- c()
for (i in 1:ncol(pairwise_combinations)){ 
  pair_tmp <- paste0(pairwise_combinations[1, i],'*',pairwise_combinations[2, i])
  pair <- c(pair,pair_tmp)
}

names(results) <- pair
saveRDS(results, "/CIMA/Result/rb_test/20250324_caQTL_rb_corr.rds")
write.csv(sample_number_df, "/CIMA/Result/rb_test/20250324_caQTL_sample_number.csv")
write.csv(sample_overlap_df, "/CIMA/Result/rb_test/20250324_caQTL_sample_overlap.csv")

library(dplyr)
library(glue)

get_independent_topsnp_pvalues <- function(a_row_of_df, phenotype_df, genotype_df, covariates_df, gene_map_dict) {
  # 提取 phenotype 和 topSNP 列
  phenotype <- a_row_of_df[2]
  phenotype <- unname(gene_map_dict[phenotype])
  topSNP_A <- a_row_of_df[3]
  topSNP_B <- a_row_of_df[4]
  
  # 从 phenotype_df 和 genotype_df 中选择所需列
  phenotype_select <- phenotype_df[, phenotype, drop = FALSE]
  genotype_select <- genotype_df[, c(topSNP_A, topSNP_B)]
  colnames(genotype_select) <- c("topsnp_CIMA", "topsnp_public")
  
  # 合并数据框
  df2 <- cbind(phenotype_select, covariates_df, genotype_select)
  df1 <- df2[, -ncol(df2)]  # 去掉最后一列
  
  # 拟合线性回归模型
  model1 <- lm(formula = formula(df1), data = df1)
  model2 <- lm(formula = formula(df2), data = df2)
  
  # 获取 model1 和 model2 的摘要
  summary_model1 <- summary(model1)
  summary_model2 <- summary(model2)
  
  # 提取 topsnp_A 的 beta 和 p 值
  beta_model1 <- summary_model1$coefficients["topsnp_CIMA", "Estimate"]
  p_value_model1 <- summary_model1$coefficients["topsnp_CIMA", "Pr(>|t|)"]
  
  beta_model2 <- summary_model2$coefficients["topsnp_CIMA", "Estimate"]
  p_value_model2 <- summary_model2$coefficients["topsnp_CIMA", "Pr(>|t|)"]
  
  # 创建一个向量存储 beta 和 p 值
  result_vector <- c(beta_model1, p_value_model1, beta_model2, p_value_model2)
  
  return(result_vector)
}

#创建命名向量，这是因为R语言会把列名的-换成.
# 假设 tss 是你的数据框
tss <- read.csv("/media/scPBMC1_AnalysisDisk1/huangzhuoli/hw5_backup/gaoyue/health/Make_Gene_TSS/tss.bed", sep = " ")
# 将 '-' 替换为 '.'，并创建新的 gene_id_map 列
tss$gene_id_map <- gsub("-", ".", tss$gene_id)
# 创建命名向量，作为一个字典的映射
gene_map_dict <- setNames( tss$gene_id_map, tss$gene_id)

#基于celltype A为中心的分析
celltypelist <- c('CD4_T','CD8_T','Myeloid','B','NK')



#处理onek1k数据
for (celltype in celltypelist) {
  print(glue('processing_{celltype}'))
  # 读取基因型
  genotype_df <- read.csv('/CIMA/Data/public_independent/genotype_sigeSNP_onek1k.csv', header = T, row.names = 1)
  
  # 读取需要检测的位点和基因
  need_to_detect <- read.csv(glue('/CIMA/Data/public_independent/20250220_onk1k_independent_input.csv'), row.names = 1)
  need_to_detect <- need_to_detect[need_to_detect[[1]] == celltype , ]
  
  # 读取表型
  phenotype_df <- read.csv(glue('/CIMA/Data/eQTL_L1/normal_dis/{celltype}.csv'), row.names = 1)
  
  # 读取协变量数据
  covariates_df <- read.csv(glue('/CIMA/Data/eQTL_L1/peer/peer_rez/factor/{celltype}.csv'), row.names = 1, header = T)
  
  # 获取样本大小
  sample_size <- nrow(covariates_df)
  
  # 确定 peer factor 列
  if (sample_size <= 150) {
    pf <- paste0("pf", 1:15)
  } else if (sample_size <= 250) {
    pf <- paste0("pf", 1:30)
  } else {
    pf <- paste0("pf", 1:35)
  }
  
  # 创建 covariates 列表
  info <- c('Age', 'sex', 'PC1', 'PC2')
  cova_need <- c(info, pf)
  
  # 从 covariates_df 中选择需要的列
  covariates_df <- covariates_df[, cova_need]
  
  # 取公共行名
  common_rownames <- Reduce(intersect, list(rownames(phenotype_df), rownames(genotype_df), rownames(covariates_df)))
  
  # 对数据框按公共行名进行匹配
  phenotype_df <- phenotype_df[match(common_rownames, rownames(phenotype_df)), ]
  genotype_df <- genotype_df[match(common_rownames, rownames(genotype_df)), ]
  covariates_df <- covariates_df[match(common_rownames, rownames(covariates_df)), ]
  
  # 应用函数处理每行数据
  apply_result <- apply(need_to_detect, 1, function(row) {
    get_independent_topsnp_pvalues(row, phenotype_df, genotype_df, covariates_df, gene_map_dict)
  })
  
  # 转置并命名列
  apply_result <- as.data.frame(t(apply_result))
  colnames(apply_result) <- c('beta_A', 'p_A', 'beta_A_regressB', 'p_A_regressB')
  
  # 合并结果
  detect_finish <- cbind(need_to_detect, apply_result)
  
  # 输出结果
  write.csv(detect_finish, file = glue('/CIMA/Result/eQTL_L1_downstream/public_independent/onek1k_{celltype}.csv'))
}

#处理immuenexut数据
for (celltype in celltypelist) {
  print(glue('processing_{celltype}'))
  # 读取基因型
  genotype_df <- read.csv('/CIMA/Data/public_independent/genotype_sigeSNP_immuenexut.csv', header = T, row.names = 1)
  
  # 读取需要检测的位点和基因
  need_to_detect <- read.csv(glue('/CIMA/Data/public_independent/20250220_immuenexut_independent_input.csv'), row.names = 1)
  need_to_detect <- need_to_detect[need_to_detect[[1]] == celltype , ]
  
  # 读取表型
  phenotype_df <- read.csv(glue('/CIMA/Data/eQTL_L1/normal_dis/{celltype}.csv'), row.names = 1)
  
  # 读取协变量数据
  covariates_df <- read.csv(glue('/CIMA/Data/eQTL_L1/peer/peer_rez/factor/{celltype}.csv'), row.names = 1, header = T)
  
  # 获取样本大小
  sample_size <- nrow(covariates_df)
  
  # 确定 peer factor 列
  if (sample_size <= 150) {
    pf <- paste0("pf", 1:15)
  } else if (sample_size <= 250) {
    pf <- paste0("pf", 1:30)
  } else {
    pf <- paste0("pf", 1:35)
  }
  
  # 创建 covariates 列表
  info <- c('Age', 'sex', 'PC1', 'PC2')
  cova_need <- c(info, pf)
  
  # 从 covariates_df 中选择需要的列
  covariates_df <- covariates_df[, cova_need]
  
  # 取公共行名
  common_rownames <- Reduce(intersect, list(rownames(phenotype_df), rownames(genotype_df), rownames(covariates_df)))
  
  # 对数据框按公共行名进行匹配
  phenotype_df <- phenotype_df[match(common_rownames, rownames(phenotype_df)), ]
  genotype_df <- genotype_df[match(common_rownames, rownames(genotype_df)), ]
  covariates_df <- covariates_df[match(common_rownames, rownames(covariates_df)), ]
  
  # 应用函数处理每行数据
  apply_result <- apply(need_to_detect, 1, function(row) {
    get_independent_topsnp_pvalues(row, phenotype_df, genotype_df, covariates_df, gene_map_dict)
  })
  
  # 转置并命名列
  apply_result <- as.data.frame(t(apply_result))
  colnames(apply_result) <- c('beta_A', 'p_A', 'beta_A_regressB', 'p_A_regressB')
  
  # 合并结果
  detect_finish <- cbind(need_to_detect, apply_result)
  
  # 输出结果
  write.csv(detect_finish, file = glue('/CIMA/Result/eQTL_L1_downstream/public_independent/immuenexut_{celltype}.csv'))
}


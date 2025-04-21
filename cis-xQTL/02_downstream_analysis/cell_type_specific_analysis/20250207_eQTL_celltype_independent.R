library(dplyr)
library(glue)

#首先定义type
type = 'eQTL'

get_independent_topsnp_pvalues <- function(a_row_of_df, phenotype_df, genotype_df, covariates_df, gene_map_dict) {
  # 提取 phenotype �? topSNP �?
  phenotype <- a_row_of_df[3]
  phenotype <- unname(gene_map_dict[phenotype])
  topSNP_A <- a_row_of_df[4]
  topSNP_B <- a_row_of_df[5]
  
  # �? phenotype_df �? genotype_df 中选择所需�?
  phenotype_select <- phenotype_df[, phenotype, drop = FALSE]
  genotype_select <- genotype_df[, c(topSNP_A, topSNP_B)]
  colnames(genotype_select) <- c("topsnp_A", "topsnp_B")
  
  # 合并数据�?
  df2 <- cbind(phenotype_select, covariates_df, genotype_select)
  df1 <- df2[, -ncol(df2)]  # 去掉最后一�?
  
  # 拟合线性回归模�?
  model1 <- lm(formula = formula(df1), data = df1)
  model2 <- lm(formula = formula(df2), data = df2)
  
  # 获取 model1 �? model2 的摘�?
  summary_model1 <- summary(model1)
  summary_model2 <- summary(model2)
  
  # 提取 topsnp_A �? beta �? p �?
  beta_model1 <- summary_model1$coefficients["topsnp_A", "Estimate"]
  p_value_model1 <- summary_model1$coefficients["topsnp_A", "Pr(>|t|)"]
  
  beta_model2 <- summary_model2$coefficients["topsnp_A", "Estimate"]
  p_value_model2 <- summary_model2$coefficients["topsnp_A", "Pr(>|t|)"]
  
  # 创建一个向量存�? beta �? p �?
  result_vector <- c(beta_model1, p_value_model1, beta_model2, p_value_model2)
  
  return(result_vector)
}

#创建命名向量，这是因为R语言会把列名�?-换成.
# 假设 tss 是你的数据框
tss <- read.csv("/CIMA/Make_Gene_TSS/tss.bed", sep = " ")
# �? '-' 替换�? '.'，并创建新的 gene_id_map �?
tss$gene_id_map <- gsub("-", ".", tss$gene_id)
# 创建命名向量，作为一个字典的映射
gene_map_dict <- setNames( tss$gene_id_map, tss$gene_id)

#基于celltype A为中心的分析
celltypelist <- read.table('/CIMA/Script/all_69_CT.txt')
celltypelist <- celltypelist$V1

for (celltype in celltypelist) {
  print(glue('processing_{celltype}'))
  # 读取基因�?
  genotype_df <- read.csv('/CIMA/Data/celltype_independent/genotype_sigeSNP.csv', header = T, row.names = 1)
  
  # 读取需要检测的位点和基�?
  need_to_detect <- read.csv(glue('/CIMA/Data/celltype_independent/{type}_df/{celltype}.csv'), row.names = 1)
  
  # 读取表型
  phenotype_df <- read.csv(glue('/CIMA/Data/{type}/normal_dis/{celltype}.csv'), row.names = 1)
  
  # 读取协变量数�?
  covariates_df <- read.csv(glue('/CIMA/Data/{type}/peer/peer_rez/factor/{celltype}.csv'), row.names = 1, header = T)
  
  # 获取样本大小
  sample_size <- nrow(covariates_df)
  
  # 确定 peer factor �?
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
  
  # �? covariates_df 中选择需要的�?
  covariates_df <- covariates_df[, cova_need]
  
  # 取公共行�?
  common_rownames <- Reduce(intersect, list(rownames(phenotype_df), rownames(genotype_df), rownames(covariates_df)))
  
  # 对数据框按公共行名进行匹�?
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
  write.csv(detect_finish, file = glue('/CIMA/Result/downstream/celltype_independent/{type}/{celltype}.csv'))
}

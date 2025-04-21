library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(glue)

eQTL_all_df = read.csv('/CIMA/Result/20250108_cis_eQTL_studywise_sig.csv')
CT_full_list = data.frame("V1" = unique(eQTL_all_df$celltype))
CT_full_list$detected_gene_num <- 0
CT_full_list$eGene_num <- 0

#统计检测的基因的数量和
detected_genelist <- list()
for (celltype in CT_full_list$V1){
  print(celltype)
  detected_df = fread(glue("/CIMA/Data/eQTL/pseudobulk/{celltype}.csv"))
  detected_df = detected_df[,-1]
  detected_number <- ncol(detected_df)
  detected_genelist[[celltype]] = colnames(detected_df)
  CT_full_list[CT_full_list$V1 == celltype,]$detected_gene_num <- detected_number
}

saveRDS(detected_genelist, file = "/CIMA/Result/summary/20250212_gene_detected_list_by_CT.rds")

eGenelist <- list()
for (celltype in CT_full_list$V1){
  print(celltype)
  eGene_df <- eQTL_all_df[eQTL_all_df$celltype == celltype,]
  eGene_number <- length(unique(eGene_df$phenotype_id))
  if (eGene_number > 0) {eGenelist[[celltype]] <- c(unique(eGene_df$phenotype_id))}
  CT_full_list[CT_full_list$V1 == celltype,]$eGene_num <- eGene_number
}

saveRDS(eGenelist, file = "/CIMA/Result/summary/20250212_eGene_list_by_CT.rds")

all_eGene <- unique(eQTL_all_df$phenotype_id)

count_gene_occurrences <- function(gene) {
  # 使用sapply检查每个基因是否在R列表的每个元素中出现，并返回逻辑向量
  gene_occurrences <- sapply(eGenelist, function(x) gene %in% x)
  # 使用sum计算逻辑向量中为TRUE的元素个数，即基因在列表中出现的次数
  count <- sum(gene_occurrences)
  #names(count) <- gene
  return(count)
}

# 使用lapply对基因列表中的每个基因应用count_gene_occurrences函数
all_eGene_cscount <- lapply(all_eGene, count_gene_occurrences)
all_eGene_cscount_df = as.data.frame.vector(all_eGene_cscount,row.names = all_eGene)
all_eGene_cscount_df$all_eGene_cscount <- as.numeric(all_eGene_cscount_df$all_eGene_cscount)
table(all_eGene_cscount_df$all_eGene_cscount == 1)
#结果是Fasle_6831_true_2769

count_gene_occurrences <- function(gene) {
  # 使用sapply检查每个基因是否在R列表的每个元素中出现，并返回逻辑向量
  gene_occurrences <- sapply(detected_genelist, function(x) gene %in% x)
  # 使用sum计算逻辑向量中为TRUE的元素个数，即基因在列表中出现的次数
  count <- sum(gene_occurrences)
  #names(count) <- gene
  return(count)
}

cs_eGene = row.names(all_eGene_cscount_df[all_eGene_cscount_df$all_eGene_cscount == 1,,drop = FALSE])
cseGene_detectedcount <- lapply(cs_eGene, count_gene_occurrences)
cseGene_detectedcount = as.data.frame.vector(cseGene_detectedcount,row.names = cs_eGene)
sum(cseGene_detectedcount == 1)
#结果是有442个只在一个celltype里面鉴定的eGene是因为只在这个celltype中被表达
  
#和泽凯的结果核对，全为true
piedf <- as.data.frame.array(table(all_eGene_cscount_df$all_eGene_cscount))
pie_df_from_zekai = read.csv("/CIMA/Result/summary/from_zekai/eGene_frequency_summary_20250117.csv",row.names = 1)
table(piedf$`table(all_eGene_cscount_df$all_eGene_cscount)` == pie_df_from_zekai$n_phenotypes)

#各个细胞类型中的cs_eGene list的数�?
cs_eGenelist <- list()
for (celltype in names(eGenelist)){
  eGenelist_cal <- eGenelist
  eGene_cs <- eGenelist_cal[[celltype]] 
  eGenelist_cal[celltype] <- NULL
  all_rest_eGene <- Reduce(union, eGenelist_cal)
  eGene_cs <- setdiff(eGene_cs, all_rest_eGene)
  cs_eGenelist[[celltype]] <- eGene_cs
}

CT_full_list$cs_eGene_num <- 0
for (celltype in names(cs_eGenelist)){
  CT_full_list[CT_full_list$V1 == celltype,"cs_eGene_num"] <- length(cs_eGenelist[[celltype]])
}

#和泽凯的结果核对,全为True
CT_full_list_zekai = read.csv('/CIMA/Result/summary/from_zekai/eGene_cellspecial_number_20250117.csv',row.names = 1)
CT_full_list$cs_eGene_num == CT_full_list_zekai$n_unique_phonetype_number
CT_full_list$eGene_num == (CT_full_list_zekai$n_unique_phonetype_number + CT_full_list_zekai$n_share_phonetype_number)


CT_full_list$cs_percentage <- CT_full_list$cs_eGene_num/CT_full_list$eGene_num
min(CT_full_list$cs_percentage)
max(CT_full_list$cs_percentage)
mean(CT_full_list$cs_percentage)

write.csv(CT_full_list,file = "/CIMA/Result/summary/20250212_eGene_number_csnumber_CT.csv")

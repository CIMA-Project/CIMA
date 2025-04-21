library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(glue)

caQTL_all_df = read.csv('/CIMA/Result/20250108_cis_caQTL_studywise_sig.csv')
CT_full_list = data.frame("V1" = unique(caQTL_all_df$celltype))
CT_full_list$detected_peak_num <- 0
CT_full_list$caPeak_num <- 0

#统计检测的基因的数量和
detected_peaklist <- list()
for (celltype in CT_full_list$V1){
  print(celltype)
  detected_df = fread(glue("/CIMA/Data/caQTL/pseudobulk/{celltype}.csv"))
  detected_df = detected_df[,-1]
  detected_number <- ncol(detected_df)
  detected_peaklist[[celltype]] = colnames(detected_df)
  CT_full_list[CT_full_list$V1 == celltype,]$detected_peak_num <- detected_number
}

saveRDS(detected_peaklist, file = "/CIMA/Result/summary/20250212_peak_detected_list_by_CT.rds")

caPeaklist <- list()
for (celltype in CT_full_list$V1){
  print(celltype)
  caPeak_df <- caQTL_all_df[caQTL_all_df$celltype == celltype,]
  caPeak_number <- length(unique(caPeak_df$phenotype_id))
  if (caPeak_number > 0) {caPeaklist[[celltype]] <- c(unique(caPeak_df$phenotype_id))}
  CT_full_list[CT_full_list$V1 == celltype,]$caPeak_num <- caPeak_number
}

saveRDS(caPeaklist, file = "/CIMA/Result/summary/20250212_caPeak_list_by_CT.rds")

all_caPeak <- unique(caQTL_all_df$phenotype_id)

count_peak_occurrences <- function(peak) {
  # 使用sapply检查每个基因是否在R列表的每个元素中出现，并返回逻辑向量
  peak_occurrences <- sapply(caPeaklist, function(x) peak %in% x)
  # 使用sum计算逻辑向量中为TRUE的元素个数，即基因在列表中出现的次数
  count <- sum(peak_occurrences)
  #names(count) <- peak
  return(count)
}

# 使用lapply对基因列表中的每个基因应用count_peak_occurrences函数
all_caPeak_cscount <- lapply(all_caPeak, count_peak_occurrences)
all_caPeak_cscount_df = as.data.frame.vector(all_caPeak_cscount,row.names = all_caPeak)
all_caPeak_cscount_df$all_caPeak_cscount <- as.numeric(all_caPeak_cscount_df$all_caPeak_cscount)
table(all_caPeak_cscount_df$all_caPeak_cscount == 1)
#结果是Fasle_23463_true_28898

count_peak_occurrences <- function(peak) {
  # 使用sapply检查每个基因是否在R列表的每个元素中出现，并返回逻辑向量
  peak_occurrences <- sapply(detected_peaklist, function(x) peak %in% x)
  # 使用sum计算逻辑向量中为TRUE的元素个数，即基因在列表中出现的次数
  count <- sum(peak_occurrences)
  #names(count) <- peak
  return(count)
}

cs_caPeak = row.names(all_caPeak_cscount_df[all_caPeak_cscount_df$all_caPeak_cscount == 1,,drop = FALSE])
cscaPeak_detectedcount <- lapply(cs_caPeak, count_peak_occurrences)
cscaPeak_detectedcount = as.data.frame.vector(cscaPeak_detectedcount,row.names = cs_caPeak)
sum(cscaPeak_detectedcount == 1)
sum(cs_caPeak %in% detected_peaklist[['cMono_CD14']])
#结果是有13027个只在一个celltype里面鉴定的caPeak是因为只在这个celltype中被表达

#和泽凯的结果核对，全为true
piedf <- as.data.frame.array(table(all_caPeak_cscount_df$all_caPeak_cscount))
pie_df_from_zekai = read.csv("/CIMA/Result/summary/from_zekai/caPeak_frequency_summary_20250117.csv",row.names = 1)
table(piedf$`table(all_caPeak_cscount_df$all_caPeak_cscount)` == pie_df_from_zekai$n_phenotypes)

#各个细胞类型中的cs_caPeak list的数�?
cs_caPeaklist <- list()
for (celltype in names(caPeaklist)){
  caPeaklist_cal <- caPeaklist
  caPeak_cs <- caPeaklist_cal[[celltype]] 
  caPeaklist_cal[celltype] <- NULL
  all_rest_caPeak <- Reduce(union, caPeaklist_cal)
  caPeak_cs <- setdiff(caPeak_cs, all_rest_caPeak)
  cs_caPeaklist[[celltype]] <- caPeak_cs
}

CT_full_list$cs_caPeak_num <- 0
for (celltype in names(cs_caPeaklist)){
  CT_full_list[CT_full_list$V1 == celltype,"cs_caPeak_num"] <- length(cs_caPeaklist[[celltype]])
}

#和泽凯的结果核对,全为True
CT_full_list_zekai = read.csv('/CIMA/Result/summary/from_zekai/caPeak_cellspecial_number_20250116.csv',row.names = 1)
CT_full_list$cs_caPeak_num == CT_full_list_zekai$n_unique_phenotypes
CT_full_list$caPeak_num == (CT_full_list_zekai$n_unique_phenotypes + CT_full_list_zekai$n_share_phenotypes)


CT_full_list$cs_percentage <- CT_full_list$cs_caPeak_num/CT_full_list$caPeak_num
min(CT_full_list$cs_percentage)
max(CT_full_list$cs_percentage)
mean(CT_full_list$cs_percentage)

write.csv(CT_full_list,file = "/CIMA/Result/summary/20250212_caPeak_number_csnumber_CT.csv")

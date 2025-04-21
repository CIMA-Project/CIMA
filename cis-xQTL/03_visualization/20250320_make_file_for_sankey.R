eGene = read.csv("/CIMA/Result/20250108_cis_eQTL_studywise_sig.csv")
eGene_share <- as.data.frame(table(eGene$phenotype_id))
eGene_detected <- readRDS("/CIMA/Result/summary/20250212_gene_detected_list_by_CT.rds")
counts <- sapply(eGene_share$Var1, function(x) sum(sapply(eGene_detected, function(y) x %in% y)))
eGene_share$Freq2 <- counts
colnames(eGene_share) <- c('gene','sig_num','detected_num')
# 使用 cut() 函数将 Freq 列分成三种类型
eGene_share$sig_num <- cut(eGene_share$sig_num ,
                     breaks = c(-Inf, 1, 10, Inf),
                     labels = c("1", "2_10", "over_10"),
                     right = TRUE)
eGene_share$detected_num <- cut(eGene_share$detected_num ,
                           breaks = c(-Inf, 1, 10, Inf),
                           labels = c("1", "2_10", "over_10"),
                           right = TRUE)
eGene_share$from_to <- paste(eGene_share$sig_num,"_xxx_",eGene_share$detected_num)
eGene_share <- as.data.frame(table(eGene_share$from_to))
split_result <- strsplit(as.character(eGene_share$Var1), "_xxx_")
eGene_share$sig_num = sapply(split_result, function(x) x[1]) # 前半部分
eGene_share$detected_num = sapply(split_result, function(x) x[2])
write.csv(eGene_share,file = '/CIMA/Result/summary/eQTL_sig_vs_detected.csv')

caPeak = read.csv("/CIMA/Result/20250108_cis_caQTL_studywise_sig.csv")
caPeak_share <- as.data.frame(table(caPeak$phenotype_id))
caPeak_detected <- readRDS("/CIMA/Result/summary/20250212_peak_detected_list_by_CT.rds")
counts <- sapply(caPeak_share$Var1, function(x) sum(sapply(caPeak_detected, function(y) x %in% y)))
caPeak_share$Freq2 <- counts
colnames(caPeak_share) <- c('gene','sig_num','detected_num')
# 使用 cut() 函数将 Freq 列分成三种类型
caPeak_share$sig_num <- cut(caPeak_share$sig_num ,
                            breaks = c(-Inf, 1, 10, Inf),
                            labels = c("1", "2_10", "over_10"),
                            right = TRUE)
caPeak_share$detected_num <- cut(caPeak_share$detected_num ,
                                 breaks = c(-Inf, 1, 10, Inf),
                                 labels = c("1", "2_10", "over_10"),
                                 right = TRUE)
caPeak_share$from_to <- paste(caPeak_share$sig_num,"_xxx_",caPeak_share$detected_num)
caPeak_share <- as.data.frame(table(caPeak_share$from_to))
split_result <- strsplit(as.character(caPeak_share$Var1), "_xxx_")
caPeak_share$sig_num = sapply(split_result, function(x) x[1]) # 前半部分
caPeak_share$detected_num = sapply(split_result, function(x) x[2])
write.csv(caPeak_share,file = '/CIMA/Result/summary/caQTL_sig_vs_detected.csv')

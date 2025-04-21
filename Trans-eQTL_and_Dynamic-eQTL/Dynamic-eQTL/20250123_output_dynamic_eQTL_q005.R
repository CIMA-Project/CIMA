library("qvalue")
library("glue")

sig_eGene_B = read.csv('/CIMA/Result/dynamic/pseudobulk/B/20250122_eGene_use_for_down_stream.csv',row.names = 1)

LRT_p_value = c()
for (n in 1:nrow(sig_eGene_B)){
  phenotype = sig_eGene_B$phenotype_id[n]
  LRT_p_add = glue('/CIMA/Result/dynamic/single_cell/B/{phenotype}.txt')
  LRT_p <- read.table(LRT_p_add)$V1[1]
  LRT_p_value = c(LRT_p_value,LRT_p)
}

sig_eGene_B$dynamic_LRT_p <- LRT_p_value
sig_eGene_B$dynamic_LRT_q <- qvalue(LRT_p_value)$qvalues
table(sig_eGene_B$dynamic_LRT_q < 0.05)
sig_eGene_B <- sig_eGene_B[sig_eGene_B$dynamic_LRT_q < 0.05,]
write.csv(sig_eGene_B,file = '/CIMA/Result/dynamic/20250123_single_cell_dynamic_QTL_Bcell.csv')
sig_eGene_B

sig_SMR_B <- read.csv('/CIMA/Result/downstream/SMR_summary/dynamic_eQTL_B.csv',row.names = 1)
sig_SMR_B <- sig_SMR_B[sig_SMR_B$topSNP %in% sig_eGene_B$variant_id,]
write.csv(sig_SMR_B,'/CIMA/Result/downstream/SMR_summary/dynamic_eQTL_B_with_sigSMR.csv')

sig_eGene_Mono= read.csv('/CIMA/Result/dynamic/pseudobulk/Monocyte/20250122_eGene_use_for_down_stream.csv',row.names = 1)

LRT_p_value = c()
for (n in 1:nrow(sig_eGene_Mono)){
  phenotype = sig_eGene_Mono$phenotype_id[n]
  LRT_p_add = glue('/CIMA/Result/dynamic/single_cell/Monocyte/{phenotype}.txt')
  LRT_p <- read.table(LRT_p_add)$V1[1]
  LRT_p_value = c(LRT_p_value,LRT_p)
}

sig_eGene_Mono$dynamic_LRT_p <- LRT_p_value
sig_eGene_Mono$dynamic_LRT_q <- qvalue(LRT_p_value)$qvalues
table(sig_eGene_Mono$dynamic_LRT_q < 0.05)
sig_eGene_Mono<- sig_eGene_Mono[sig_eGene_Mono$dynamic_LRT_q < 0.05,]
write.csv(sig_eGene_Mono,file = '/CIMA/Result/dynamic/20250123_single_cell_dynamic_QTL_Monocyte.csv')

sig_SMR_Mono <- read.csv('/CIMA/Result/downstream/SMR_summary/dynamic_eQTL_Monocyte.csv',row.names = 1)
sig_SMR_Mono  <- sig_SMR_Mono[sig_SMR_Mono$topSNP %in% sig_eGene_Mono$variant_id,]
write.csv(sig_SMR_Mono,'/CIMA/Result/downstream/SMR_summary/dynamic_eQTL_Mono_with_sigSMR.csv')

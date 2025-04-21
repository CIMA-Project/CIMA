set.seed(8888)
setwd('/CIMA/Data/GWAS/')
library(tidyverse)
library(data.table)

#输出代谢物的数据

metabolity <- fread('./Metabolites.csv') %>%
  `colnames<-`( gsub('-1', '', colnames(.)) )
sum(is.na(metabolity))

geno.sample <- fread('./GWAS_lip_meta_378.txt',header = F)
geno.sample$V2 <- geno.sample$V1
colnames(geno.sample) <- c('FID','IID')

metabolity.name <- metabolity$SAMPLE[-1]

# phenotype files
meta.list <- matrix(NA, nrow = length(metabolity.name), ncol = 2)
for (meta.single in metabolity.name) {
  
  num <- match(meta.single, metabolity.name)
  meta.df <- metabolity %>%
    filter(SAMPLE == meta.single) %>%
    t() %>% as.data.frame() %>%
    rownames_to_column('sample')
  
  meta.df <- meta.df[-1, ]
  meta.df[ ,2] <- as.numeric(meta.df[ ,2])
  
  print(head(meta.df))
  
  df <- geno.sample%>%
    inner_join(meta.df, by = c('FID' = 'sample'))
  
  df <- df %>%
    mutate( norm = qnorm( (rank(.$V1, na.last = 'keep') - 0.5)/sum(!is.na(.$V1)) ) ) %>%
    dplyr::select(-V1)
  
  fwrite(df, sprintf('/CIMA/Data/GWAS/meta/META_%s.txt', num), sep = '\t', row.names = F, col.names = F)
  meta.list[num, ] <- c(meta.single, sprintf('META_%s', num))
  
}

#验证和之前输入的相似程度
test_df = read.table(sprintf(/CIMA/GWAS/Metabolites/pheno/META_%s.txt', num))
sum(df$FID %in% test_df$V1)
test_df_sorted <- test_df[match(df$FID, test_df$V1), ]
print(cor(df$norm,test_df_sorted$V3))

meta.list %>%
  fwrite('./meta.list.csv', col.names = F)

#生成qcovar(pc1,pc2和年龄)和covar(性别和批次)
sample_info = read.csv('/CIMA/Data/sample_info_new.csv')
sample_info = sample_info[match(geno.sample$FID,sample_info$BGE_name),]
sum(sample_info$BGE_name == geno.sample$FID)
qcovar = sample_info[,c('BGE_name','BGE_name','PC1','PC2','Age')]
covar = sample_info[,c('BGE_name','BGE_name','sex')]
batch <- unlist(metabolity[1,.SD])
names(batch) <- colnames(metabolity)
covar$batch <- batch[covar$BGE_name]

#测试和之前的相同程度
covar_test= read.table(/CIMA/GWAS/Metabolites/covar.3.txt')
covar_test=covar_test[match(covar$BGE_name,covar_test$V1),]
sum(covar_test$V3 == covar$sex)
sum(covar_test$V4 == covar$batch)
qcovar_test = read.table(/CIMA/GWAS/Metabolites/qcovar.txt')
qcovar_test=qcovar_test[match(covar$BGE_name,qcovar_test$V1),]
print(cor(qcovar$PC1,qcovar_test$V3))
print(cor(qcovar$PC2,qcovar_test$V4))
sum(qcovar$Age == qcovar_test$V5)
#测试结果是完全相同

#输出协变量文件
fwrite(qcovar, './qcovar_meta_lipid.txt', sep = '\t', col.names = F)
fwrite(covar, './covar_meta.txt', sep = '\t', col.names = F)
fwrite(covar[1:3], './covar_lipid.txt', sep = '\t', col.names = F)


##处理脂质组
lipid = fread('./lipidome.tsv')
#lipid = lipid %>%
#  filter(sample %in% geno.sample$FID)
lipid = lipid[match(geno.sample$FID,lipid$sample),]
sum(lipid$sample == geno.sample$FID)

na.sum = matrix(NA, nrow = ncol(lipid), ncol = 2)
for (i in 1:ncol(lipid)) {
  na.sum[i,] = c(names(lipid)[i], sum(is.na(as.data.frame(lipid)[ ,i])))
  # print(names(lipid)[i])
  # print(sum(is.na(as.data.frame(lipid)[ ,i])))
}

hist(as.numeric(na.sum[,2]), breaks = seq(0,500,5), freq = T, xlab = 'No. of NA value', main = NULL)

keep.lipid = na.sum[,1][as.numeric(na.sum[,2]) <= 378*0.25]

#测试一下上一版本的显著结果是否在里面
c('LPC 20:4','LPE 22:6','PE 18:2-20:1','TAG48:4-FA12:0','TAG48:4-FA18:2') %in% keep.lipid 

lipid = as.data.frame(lipid)[,c(keep.lipid)]
zzz = names(lipid)[-1]
lipid.list <- matrix(NA, nrow = length(zzz), ncol = 2)
for (i in 1:length(zzz)) {
  
  lipid.zzz = lipid %>%
    dplyr::select(sample, all_of(zzz[i]))
  colnames(lipid.zzz)[2] = 'zzz'
  
  ## half minimize imputation
  min.zzz = min(lipid.zzz$zzz, na.rm = T)
  lipid.zzz[is.na(lipid.zzz)] = (min.zzz / 2)
  
  ## rank inverse normalization
  lipid.norm = lipid.zzz %>%
    mutate(norm = qnorm( (rank(.$zzz, na.last = 'keep') - 0.5) / sum(!is.na(.$zzz)) )) %>%
    mutate(FID = .$sample) %>%
    dplyr::select(FID, sample, norm)
  
  #norm.test.mat[i,] = c(zzz[i], shapiro.test(lipid.norm$norm)$p.value)
  
  fwrite(lipid.norm, sprintf('./lipid/lipid_%s.txt', i), col.names = F, sep = '\t')
  lipid.list[i, ] <- c(zzz[i], sprintf('lipid_%s', i))
}

lipid.list %>%
  fwrite('./lipid.list.csv', col.names = F)

#验证和之前输入的相似程度
test_df = read.table(/CIMA/GWAS/Lipid/pheno/lipid_1.txt')
test_df2 = read.table('./lipid/lipid_1.txt')
test_df_sorted <- test_df[match(test_df2$V1, test_df$V1), ]
print(cor(test_df2$V3,test_df_sorted$V3))

meta.list %>%
  fwrite('./meta.list.csv', col.names = F)
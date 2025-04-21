library(lme4)
library(Matrix)
library(glue)

args <- commandArgs(trailingOnly = TRUE)
celltype_gene <- args[1]
print(celltype_gene)
celltype <- strsplit(celltype_gene, "_xxx_")[[1]][1]
gene <- strsplit(celltype_gene, "_xxx_")[[1]][2]

phenotype <- read.csv(glue("/CIMA/Data/dynamic/single_cell/{celltype}/{gene}.csv"),row.names = 1)
cov_df <- read.csv(glue("/CIMA/Data/dynamic/single_cell/{celltype}_cov_df.csv"),row.names = 1)
data = cbind(phenotype,cov_df)
data$total_counts <- scale(log(data$total_counts))

full_model <- lme4::glmer(formula = E ~ G + (1 | sample) + Age + sex + total_counts +pct_counts_mt + PC1 + PC2  + mRNA_PC1 + mRNA_PC2 + mRNA_PC3 + mRNA_PC4 + mRNA_PC5 + ptime + G*ptime, 
                          family = "poisson", nAGQ = 0, data= data, control = glmerControl(optimizer = "nloptwrap"))

null_model <- lme4::glmer(formula = E ~ G + (1 | sample) + Age + sex + total_counts +pct_counts_mt + PC1 + PC2  + mRNA_PC1 + mRNA_PC2 + mRNA_PC3 + mRNA_PC4 + mRNA_PC5 + ptime, 
                          family = "poisson", nAGQ = 0, data= data, control = glmerControl(optimizer = "nloptwrap"))

model_lrt <- anova(null_model, full_model)
write(model_lrt$`Pr(>Chisq)`[2], file = glue("/CIMA/Result/dynamic/single_cell/{celltype}/{gene}.txt"))
print(glue("{celltype_gene}_done"))
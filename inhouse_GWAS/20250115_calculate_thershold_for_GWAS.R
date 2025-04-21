set.seed(8888)
setwd('/CIMA/Data/GWAS/')
library(tidyverse)
library(data.table)

#检测代谢物之间的相关性
metabolite.all <- list()
for (i in 1:321) {
  
  metabolite.all[[i]] <- fread(sprintf('./meta/META_%s.txt', i)) %>%
    dplyr::select(-V2) %>%
    `rownames<-`(.$V1) %>%
    dplyr::select(-V1) %>%
    `colnames<-`(c(paste0('meta', i)))
  
}
metabolite.mat <- bind_cols(metabolite.all)
pca <- prcomp(x = metabolite.mat, center = F, scale. = F)
explained_variance_ratio <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
which.max(explained_variance_ratio >= 0.95)
explained_variance_ratio[148]
meta_thershold = 5e-8/148
#阈值是3.38e-10

#鉴定独立的脂质
lipid.all <- list()
for (i in 1:653) {
  
  lipid.all[[i]] <- fread(sprintf('./lipid/lipid_%s.txt', i)) %>%
    dplyr::select(-V2) %>%
    `rownames<-`(.$V1) %>%
    dplyr::select(-V1) %>%
    #         `colnames<-`(c('sample', paste0('lipid', i)))
    `colnames<-`(c(paste0('lipid', i)))
  
}
lipid.mat <- bind_cols(lipid.all)
pca <- prcomp(x = lipid.mat, center = F, scale. = F)
explained_variance_ratio <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
which.max(explained_variance_ratio >= 0.95)

cov_m = cov(lipid.mat)
ev = eigen(cov_m)
ev_val = ev$values
sum(ev_val)^2
sum(ev_val^2)
sum(ev_val)^2 / sum(ev_val^2)

#sum(ev_val)^2 / sum(ev_val^2) 的值是6.3,近似等于6
lipid_thershold = 5e-8/6

#检测细胞比例之间的相关性
#检测代谢物之间的相关性
celltype.all <- list()
for (i in 1:73) {
  
  celltype.all[[i]] <- fread(sprintf('./cell_propotion/celltype_%s.txt', i)) %>%
    dplyr::select(-V2) %>%
    `rownames<-`(.$V1) %>%
    dplyr::select(-V1) %>%
    `colnames<-`(c(paste0('celltype', i)))
  
}

celltype.mat <- bind_cols(celltype.all)
pca <- prcomp(x = celltype.mat, center = F, scale. = F)
explained_variance_ratio <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
which.max(explained_variance_ratio >= 0.95)
explained_variance_ratio[47]
celltype_thershold = 5e-8/47
5e-8/73

cov_m = cov(celltype.mat)
ev = eigen(cov_m)
ev_val = ev$values
sum(ev_val)^2
sum(ev_val^2)
sum(ev_val)^2 / sum(ev_val^2)
5e-8/20

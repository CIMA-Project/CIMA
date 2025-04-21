##################################################################################_
# 241231
# huangzhuoli 
# Health Project
# PEER Factor
# 2000 HVG
# Environment: peer

##################################################################################_
setwd('/CIMA/Data/dynamic/pseudobulk')

library(peer)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
cell <- args[1]
print(cell)

expr <- fread(sprintf('./top2000_normal_dis/%s.csv', cell))
cova_mat <- fread('../../sample_info_new.csv')
cova_mat <- cova_mat[cova_mat$BGE_name %in% expr$V1,]
cova_mat <- cova_mat[match(expr$V1, cova_mat$BGE_name), ]

expr_mat <- expr[,-1]
cova_mat <- cova_mat[,-1]
##### PEER-------------------------------------------------------------------------
# Set PEER paramaters based on the instructions from PEER package website
model = PEER()
PEER_setPhenoMean(model, as.matrix(expr_mat))
PEER_setCovariates(model, as.matrix(cova_mat))
dim(PEER_getPhenoMean(model))
# Set to generate 10 PEER factors: K=10
K <- 50
PEER_setNk(model,K)
PEER_getNk(model)
# Calculate
PEER_update(model)
##### Observing output--------------------------------------------------------------
pf.names <- paste('pf',1:K, sep = '')
# factors: NxK matrix
f <- PEER_getX(model)
f_df <- data.frame(f)
colnames(f_df) <- c(colnames(cova_mat),pf.names)
f_df$sample <- expr$V1
f_df <- f_df[,c(ncol(f_df),1:(ncol(f_df)-1))]
fwrite(f_df, sprintf('./peer/peer_rez/factor/%s.csv',cell))

# weights: GxK matrix
w <- PEER_getW(model)
w_df <- data.frame(w)
colnames(w_df) <- c(colnames(cova_mat),pf.names)
w_df$geneid <- colnames(expr_mat)
w_df <- w_df[,c(ncol(w_df),1:(ncol(w_df)-1))]
fwrite(w_df, sprintf('./peer/peer_rez/weights/%s.csv',cell))

# precision: Kx1 matrix
p <- PEER_getAlpha(model)
p_df <- data.frame(p)
p_df$peer_factor <- c(colnames(cova_mat),pf.names)
p_df <- p_df[,c(2,1)]
fwrite(p_df, sprintf('./peer/peer_rez/precision/%s.csv',cell))

# residuals: NxG matrix
r <- PEER_getResiduals(model)
r_df <- data.frame(r)
colnames(r_df) <- colnames(expr_mat)
r_df <- cbind(sample = expr$V1, r_df)
fwrite(r_df, sprintf('./peer/peer_rez/residuals/%s.csv',cell))

print(sprintf('%s DONE', cell))

quit()



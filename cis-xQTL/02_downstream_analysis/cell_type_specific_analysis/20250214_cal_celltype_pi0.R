library(qvalue)
library(dplyr)
library(ggplot2)
library(data.table)
library(glue)

#对qvalue包的pi0est函数进行修改,已经在能用qvalue包跑出的结果进行验证结果一�?
pi0est1 <- function(p, lambda = seq(0.05,0.95,0.05), pi0.method = "smoother",
                    smooth.df = 3, smooth.log.pi0 = FALSE, ...) {
  ll <- length(lambda)
  if (ll == 1) {
    pi0 <- mean(p >= lambda)/(1 - lambda)
    pi0.lambda <- pi0
    pi0 <- min(pi0, 1)
    pi0Smooth <- NULL
  } else {
    ind <- length(lambda):1
    count<-tabulate(findInterval(p, vec=lambda))[ind]
    # 将元素（如果�? NA）替换为 0,将第一个元素（如果�? NA）替换为 0,bug出现在此�?,tabulate()函数的限制导致，可参考https://github.com/StoreyLab/qvalue/issues/38
    count[is.na(count)] <- 0
    pi0 <- cumsum(count) / (length(p) * (1-lambda[ind]))
    pi0 <- pi0[ind]
    pi0.lambda <- pi0
    # Smoother method approximation
    if (pi0.method == "smoother") {
      if (smooth.log.pi0) {
        pi0 <- log(pi0)
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- exp(predict(spi0, x = lambda)$y)
        pi0 <- min(pi0Smooth[ll], 1)
      } else {
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- predict(spi0, x = lambda)$y
        pi0 <- min(pi0Smooth[ll], 1)
      }
    }
  }
  return(pi0)
}

##产生eQTL__pi1数据
DT_eQTL = fread('/CIMA/Result/pi1_test/20250213_celltype_pi1_eQTLdata.csv')
result <- DT_eQTL[, pi0est1(pval_nominal), by = celltype_pair]
result <- result[, c("celltype.s", "celltype.compared") := tstrsplit(celltype_pair, "_xxx_", fixed = TRUE)]
write.csv(result,'/CIMA/Result/pi1_test/celltype/20250214_eGene_pi0.csv')

result_compare = read.csv('/CIMA/Result/pi1_test/from_zekai/final_eGene_pi0_results_20250122.csv')
result_compare$celltype_pair <- paste0(result_compare$celltype.s, "_xxx_", result_compare$celltype.compared)
# 合并两个数据框，�? celltype_pair �?
merged_df <- merge(result, result_compare, by = "celltype_pair", suffixes = c("_df1", "_df2"))
sum(round(merged_df$V1,3) == round(merged_df$pi0,3))


##产生caQTL__pi1数据
DT_caQTL = fread('/CIMA/Result/pi1_test/20250213_celltype_pi1_caQTLdata.csv')
result <- DT_caQTL[, pi0est1(pval_nominal), by = celltype_pair]
result <- result[, c("celltype.s", "celltype.compared") := tstrsplit(celltype_pair, "_xxx_", fixed = TRUE)]
write.csv(result,'/CIMA/Result/pi1_test/celltype/20250214_caPeak_pi0.csv')

result_compare = read.csv('/CIMA/Result/pi1_test/from_zekai/final_caPeak_pi0_results_20250126.csv')
result_compare$celltype_pair <- paste0(result_compare$celltype.s, "_xxx_", result_compare$celltype.compared)
# 合并两个数据框，�? celltype_pair �?
merged_df <- merge(result, result_compare, by = "celltype_pair", suffixes = c("_df1", "_df2"))
sum(round(merged_df$V1,3) == round(merged_df$pi0,3))

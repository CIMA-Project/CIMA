library(qvalue)
library(dplyr)
library(ggplot2)
library(data.table)
library(glue)

#对qvalue包的pi0est函数进行修改,已经在能用qvalue包跑出的结果进行验证结果一致
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
    # 将元素（如果有 NA）替换为 0,将第一个元素（如果有 NA）替换为 0,bug出现在此处,tabulate()函数的限制导致，可参考https://github.com/StoreyLab/qvalue/issues/38
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

##产生onek1k__pi0数据
DT_eQTL = fread('/CIMA/Result/eQTL_L1_downstream/20250219_onek1k_vs_own.csv')
result <- DT_eQTL[, pi0est1(pval_nominal), by = celltype_pair]
result <- result[, c("celltype.s", "celltype.compared") := tstrsplit(celltype_pair, "_xxx_", fixed = TRUE)]
write.csv(result,'/CIMA/Result/pi1_test/public/20250219_onek1k_pi0.csv')


##产生immuenexut__pi1数据
DT_caQTL = fread('/CIMA/Result/eQTL_L1_downstream/20250219_immuenexut_vs_own.csv')
result <- DT_caQTL[, pi0est1(pval_nominal), by = celltype_pair]
result <- result[, c("celltype.s", "celltype.compared") := tstrsplit(celltype_pair, "_xxx_", fixed = TRUE)]
write.csv(result,'/CIMA/Result/pi1_test/public/20250219_immuenexut_pi0.csv')


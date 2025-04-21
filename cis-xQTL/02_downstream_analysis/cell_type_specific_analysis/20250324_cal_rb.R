library("data.table")
# calcualte rb
# b1 and se1 represent the estimate and SE for eQTLs across probes in one tissue, b2 and se2 represent the estimate and SE for eQTLs across probes in the other tissue
# theta = sample overlap * phnotypic correlation; theta could also be estimated from null SNPs; if two samples were independent, theta = 0.
# Please note that the effect allele of one SNP between two tissues should be the same.
calcu_cor_true<-function(b1,se1,b2,se2,theta){
  
  idx=which(is.infinite(b1) | is.infinite(b2) | is.infinite(se1) | is.infinite(se2));
  if(length(idx)>0){
    b1=b1[-idx];se1=se1[-idx]
    b2=b2[-idx];se2=se2[-idx]
    theta=theta[-idx]
  }
  
  var_b1=var(b1,na.rm=T)-mean(se1^2,na.rm=T)
  var_b2=var(b2,na.rm=T)-mean(se2^2,na.rm=T)
  if(var_b1<0){
    var_b1=var(b1,na.rm=T)
  }
  if(var_b2<0){
    var_b2=var(b2,na.rm=T)
  }
  cov_b1_b2=cov(b1,b2,use="complete.obs")-mean(theta,na.rm=T)*sqrt(mean(se1^2,na.rm=T)*mean(se2^2,na.rm=T))
  r=cov_b1_b2/sqrt(var_b1*var_b2)
  
  r_jack=c()
  n=length(b1)
  for(k in 1:n) {
    b1_jack=b1[-k];se1_jack=se1[-k];var_b1_jack=var(b1_jack,na.rm=T)-mean(se1_jack^2,na.rm=T)
    b2_jack=b2[-k];se2_jack=se2[-k];var_b2_jack=var(b2_jack,na.rm=T)-mean(se2_jack^2,na.rm=T)
    if(var_b1_jack<0){
      var_b1_jack=var(b1_jack,na.rm=T);
    }
    if(var_b2_jack<0){
      var_b2_jack=var(b2_jack,na.rm=T);
    }
    theta_jack=theta[-k];
    cov_e1_jack_e2_jack=mean(theta_jack,na.rm=T)*sqrt(mean(se1_jack^2,na.rm=T)*mean(se2_jack^2,na.rm=T))
    cov_b1_b2_jack=cov(b1_jack,b2_jack,use="complete.obs")-cov_e1_jack_e2_jack
    r_tmp=cov_b1_b2_jack/sqrt(var_b1_jack*var_b2_jack)
    r_jack=c(r_jack,r_tmp)
  }
  r_mean=mean(r_jack,na.rm=T)
  idx=which(is.na(r_jack))
  if(length(idx)>0){
    se_r=sqrt((n-1)/n*sum((r_jack[-idx]-r_mean)^2))
  }else{
    se_r=sqrt((n-1)/n*sum((r_jack-r_mean)^2))
  }
  res<-cbind(r,se_r)
  return(res)
}


theta_df = read.csv('/CIMA/Result/rb_test/20250324_eQTL_sample_overlap.csv',row.names = 1,header = T)
theta_corr_list = readRDS('/CIMA/Result/rb_test/20250324_eQTL_rb_corr.rds')

rb_df = theta_df
rb_df[,] = NA
rb_se_df = rb_df

beta_se_df = fread('/CIMA/Result/rb_test/20250324_eQTL_beta_se.csv')
beta_se_df[, phenotype_id := tstrsplit(pair, "_", keep = 1)]

for (CT_ref in row.names(rb_df)){
  for (CT_query in row.names(rb_df)[row.names(rb_df) != CT_ref]){
    
    celltype_pair_use <- paste0(CT_ref,"_xxx_",CT_query)
    beta_se_df_select  <- beta_se_df [celltype_pair == celltype_pair_use]
    
    matching <- c(paste0(CT_ref,"*",CT_query),paste0(CT_query,"*",CT_ref))[c(paste0(CT_ref,"*",CT_query),paste0(CT_query,"*",CT_ref)) %in% names(theta_corr_list)]
    
    if (theta_df[CT_ref,gsub('-','.',CT_query)] == 1){theta_final <- 1} 
    else {
      interect_gene = intersect(gsub('-','.',beta_se_df_select[,phenotype_id]),names(theta_corr_list[[matching]]))
      theta_corr <- theta_corr_list[[matching]][interect_gene]
      theta_final <- theta_df[CT_ref,gsub('-','.',CT_query)]*theta_corr
    }
    
    rb_result <- calcu_cor_true(b1 = beta_se_df_select[,slope_A],se1 = beta_se_df_select[,slope_se_A],b2 = beta_se_df_select[,slope_B],se2 =  beta_se_df_select[,slope_se_B],theta = theta_final)
    
    rb_df[match(CT_ref,row.names(rb_df)),gsub('-','.',CT_query)] <- rb_result[1,1]
    rb_se_df[match(CT_ref,row.names(rb_df)),gsub('-','.',CT_query)] <- rb_result[1,2]
  }
}

write.csv(rb_df,file = '/CIMA/Result/rb_test/20250324_eQTL_rb.csv')
write.csv(rb_se_df,file = '/CIMA/Result/rb_test/20250324_eQTL_rb_se.csv')


#caQTL
theta_df = read.csv('/CIMA/Result/rb_test/20250324_caQTL_sample_overlap.csv',row.names = 1,header = T)
theta_corr_list = readRDS('/CIMA/Result/rb_test/20250324_caQTL_rb_corr.rds')

rb_df = theta_df
rb_df[,] = NA
rb_se_df = rb_df

beta_se_df = fread('/CIMA/Result/rb_test/20250324_caQTL_beta_se.csv')
beta_se_df[, phenotype_id := tstrsplit(pair, "_", keep = 1)]

for (CT_ref in row.names(rb_df)){
  for (CT_query in row.names(rb_df)[row.names(rb_df) != CT_ref]){
    
    celltype_pair_use <- paste0(CT_ref,"_xxx_",CT_query)
    beta_se_df_select  <- beta_se_df [celltype_pair == celltype_pair_use]
    
    matching <- c(paste0(CT_ref,"*",CT_query),paste0(CT_query,"*",CT_ref))[c(paste0(CT_ref,"*",CT_query),paste0(CT_query,"*",CT_ref)) %in% names(theta_corr_list)]
    
    if (theta_df[CT_ref,gsub('-','.',CT_query)] == 1){theta_final <- 1} else {
      interect_gene = intersect(gsub('[:|-]','.',beta_se_df_select[,phenotype_id]),names(theta_corr_list[[matching]]))
      theta_corr <- theta_corr_list[[matching]][interect_gene]
      theta_final <- theta_df[CT_ref,gsub('-','.',CT_query)]*theta_corr
    }
    
    rb_result <- calcu_cor_true(b1 = beta_se_df_select[,slope_A],se1 = beta_se_df_select[,slope_se_A],b2 = beta_se_df_select[,slope_B],se2 =  beta_se_df_select[,slope_se_B],theta = theta_final)
    
    rb_df[match(CT_ref,row.names(rb_df)),gsub('-','.',CT_query)] <- rb_result[1,1]
    rb_se_df[match(CT_ref,row.names(rb_df)),gsub('-','.',CT_query)] <- rb_result[1,2]
  }
}

write.csv(rb_df,file = '/CIMA/Result/rb_test/20250324_caQTL_rb.csv')
write.csv(rb_se_df,file = '/CIMA/Result/rb_test/20250324_caQTL_rb_se.csv')

min_index <- which.min(as.matrix(rb_df))
min_position <- arrayInd(min_index, dim(rb_df))
row.names(rb_df)[min_position[1]]
colnames(rb_df)[min_position[2]]

                 
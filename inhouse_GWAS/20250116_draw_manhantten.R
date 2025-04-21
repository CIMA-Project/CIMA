library("CMplot")
library('data.table')
library("dplyr")
setwd('/CIMA/Result/plot/')
lipid = read.csv('/CIMA/Result/GWAS/sig_lipid.csv')
lipid_id = lipid$V2
metabolism = read.csv('/CIMA/Result/GWAS/sig_meta.csv')
metabolism_id = metabolism$V2

gwas_df = fread(paste0('/CIMA/CIMA_r1/Result/GWAS_lipid/',lipid_id[1],'_378.mlma'),sep = '\t')
gwas_df = gwas_df[,c('SNP','Chr','bp')]
colnames(gwas_df) = c('SNP','Chromosome','Position')

for (each_lipid in lipid_id) {
  temp_GWAS_df = fread(paste0('/CIMA/CIMA_r1/Result/GWAS_lipid/',each_lipid,'_378.mlma'),sep = '\t')
  gwas_df[,lipid[lipid$V2 == each_lipid,"V1"]] = temp_GWAS_df$p
}

sig_df <- gwas_df %>%
  filter(if_any(4:7, ~ . < 8.3e-9))
no_sig_df <- gwas_df %>%
  filter(if_all(4:7, ~ . >= 8.3e-9))

sampled_df1 <- no_sig_df %>% sample_frac(0.02)
df_plot <- rbind(sig_df,sampled_df1)

setwd('/CIMA/Result/plot/')
CMplot(df_plot, plot.type="m",col=c("gray50"),multraits=TRUE,threshold=8.3e-9,threshold.lty=1,
        threshold.lwd=c(1,1), threshold.col=c("black","grey"),amplify=TRUE,
        chr.den.col=NULL, signal.col=c("#011627","#F71735","#41EAD4","#FF9F1C"),signal.cex=1, 
        file="jpg",file.name="20250116_lipid",dpi=300,file.output=TRUE,verbose=TRUE,
        points.alpha=225,legend.ncol=3, legend.pos="middle")

CMplot(df_plot, plot.type="m",col=c("gray50"),multraits=TRUE,threshold=8.3e-9,threshold.lty=1,
       threshold.lwd=c(1,1), threshold.col=c("black","grey"),amplify=TRUE,
       chr.den.col=NULL, signal.col=c("#011627","#F71735","#41EAD4","#FF9F1C"),signal.cex=1, 
       file="pdf",file.name="20250116_lipid",dpi=300,file.output=TRUE,verbose=TRUE,
       points.alpha=225,legend.ncol=3, legend.pos="middle")

CMplot(df_plot,plot.type="q",col=c("#011627","#F71735","#41EAD4","#FF9F1C"),multraits=TRUE,ylab.pos=2,signal.pch=c(19,6,4),signal.cex=1.2,signal.col="red",
       conf.int=TRUE,box=FALSE,axis.cex=1,file="jpg",file.name="20250116_lipid_qq",dpi=300,file.output=TRUE,
       verbose=TRUE,width=8,height=8)

CMplot(df_plot,plot.type="q",col=c("#011627","#F71735","#41EAD4","#FF9F1C"),multraits=TRUE,ylab.pos=2,signal.pch=c(19,6,4),signal.cex=1.2,signal.col="red",
       conf.int=TRUE,box=FALSE,axis.cex=1,file="pdf",file.name="20250116_lipid_qq",dpi=300,file.output=TRUE,
       verbose=TRUE,width=8,height=8)

write.csv(sig_df,file = '/CIMA/Result/GWAS/20250116_lipid_sig_snp.csv')

##########################for legend############################################

sampled_df1 <- no_sig_df %>% sample_frac(0.00005)
df_plot <- rbind(sig_df,sampled_df1)

setwd('/CIMA/Result/plot/')
CMplot(df_plot, plot.type="m",col=c("gray50"),multraits=TRUE,threshold=8.3e-9,threshold.lty=1,
       threshold.lwd=c(1,1), threshold.col=c("black","grey"),amplify=TRUE,
       chr.den.col=NULL, signal.col=c("#011627","#F71735","#41EAD4","#FF9F1C"),signal.cex=1, 
       file="jpg",file.name="20250116_lipid_for_legend",dpi=300,file.output=TRUE,verbose=TRUE,
       points.alpha=225,legend.ncol=3, legend.pos="middle")

CMplot(df_plot, plot.type="m",col=c("gray50"),multraits=TRUE,threshold=8.3e-9,threshold.lty=1,
       threshold.lwd=c(1,1), threshold.col=c("black","grey"),amplify=TRUE,
       chr.den.col=NULL, signal.col=c("#011627","#F71735","#41EAD4","#FF9F1C"),signal.cex=1, 
       file="pdf",file.name="20250116_lipid_for_legend",dpi=300,file.output=TRUE,verbose=TRUE,
       points.alpha=225,legend.ncol=3, legend.pos="middle")

CMplot(df_plot,plot.type="q",col=c("#011627","#F71735","#41EAD4","#FF9F1C"),multraits=TRUE,ylab.pos=2,signal.pch=c(19,6,4),signal.cex=1.2,signal.col="red",
       conf.int=TRUE,box=FALSE,axis.cex=1,file="jpg",file.name="20250116_lipid_qq_for_legend",dpi=300,file.output=TRUE,
       verbose=TRUE,width=8,height=8)

CMplot(df_plot,plot.type="q",col=c("#011627","#F71735","#41EAD4","#FF9F1C"),multraits=TRUE,ylab.pos=2,signal.pch=c(19,6,4),signal.cex=1.2,signal.col="red",
       conf.int=TRUE,box=FALSE,axis.cex=1,file="pdf",file.name="20250116_lipid_qq_for_legend",dpi=300,file.output=TRUE,
       verbose=TRUE,width=8,height=8)


##########################metabolism############################################
gwas_df = fread(paste0('/CIMA/CIMA_r1/Result/GWAS_meta/',metabolism_id[1],'_378.mlma'),sep = '\t')
gwas_df = gwas_df[,c('SNP','Chr','bp')]
colnames(gwas_df) = c('SNP','Chromosome','Position')

for (each_metabo in metabolism_id) {
  temp_GWAS_df = fread(paste0('/CIMA/CIMA_r1/Result/GWAS_meta/',each_metabo,'_378.mlma'),sep = '\t')
  gwas_df[,metabolism[metabolism$V2 == each_metabo,"V1"]] = temp_GWAS_df$p
}

5e-8/148

sig_df <- gwas_df %>%
  filter(if_any(4:10, ~ . < 3.378378e-10))
no_sig_df <- gwas_df %>%
  filter(if_all(4:10, ~ . >= 3.378378e-10))

sampled_df1 <- no_sig_df %>% sample_frac(0.02)
df_plot <- rbind(sig_df,sampled_df1)

CMplot(df_plot, plot.type="m",col=c("gray50"),multraits=TRUE,threshold=3.378378e-10,threshold.lty=1,
       threshold.lwd=c(1,1), threshold.col=c("black","grey"),amplify=TRUE,
       chr.den.col=NULL, signal.col=c("#011627","#F71735","#41EAD4","#6883BA","#FF9F1C","#118AB2","#06D6A0"),signal.cex=1, 
       file="jpg",file.name="20250116_metabolism",dpi=300,file.output=TRUE,verbose=TRUE,
       points.alpha=225,legend.ncol=3, legend.pos="middle")

CMplot(df_plot, plot.type="m",col=c("gray50"),multraits=TRUE,threshold=3.378378e-10,threshold.lty=1,
       threshold.lwd=c(1,1), threshold.col=c("black","grey"),amplify=TRUE,
       chr.den.col=NULL, signal.col=c("#011627","#F71735","#41EAD4","#6883BA","#FF9F1C","#118AB2","#06D6A0"),signal.cex=1, 
       file="pdf",file.name="20250116_metabolism",dpi=300,file.output=TRUE,verbose=TRUE,
       points.alpha=225,legend.ncol=3, legend.pos="middle")

CMplot(df_plot,plot.type="q",col=c("#011627","#F71735","#41EAD4","#6883BA","#FF9F1C","#118AB2","#06D6A0"),multraits=TRUE,ylab.pos=2,signal.pch=c(19,6,4),signal.cex=1.2,signal.col="red",
         conf.int=TRUE,box=FALSE,axis.cex=1,file="jpg",file.name="20250116_metabolism_qq",dpi=300,file.output=TRUE,
         verbose=TRUE,width=8,height=8)

CMplot(df_plot,plot.type="q",col=c("#011627","#F71735","#41EAD4","#6883BA","#FF9F1C","#118AB2","#06D6A0"),multraits=TRUE,ylab.pos=2,signal.pch=c(19,6,4),signal.cex=1.2,signal.col="red",
       conf.int=TRUE,box=FALSE,axis.cex=1,file="pdf",file.name="20250116_metabolism_qq",dpi=300,file.output=TRUE,
       verbose=TRUE,width=8,height=8)


write.csv(sig_df,file = '/CIMA/Result/GWAS/20250116_meta_sig_snp.csv')
##########################for legend############################################

sampled_df1 <- no_sig_df %>% sample_frac(0.00005)
df_plot <- rbind(sig_df,sampled_df1)

CMplot(df_plot, plot.type="m",col=c("gray50"),multraits=TRUE,threshold=3.378378e-10,threshold.lty=1,
       threshold.lwd=c(1,1), threshold.col=c("black","grey"),amplify=TRUE,
       chr.den.col=NULL, signal.col=c("#011627","#F71735","#41EAD4","#6883BA","#FF9F1C","#118AB2","#06D6A0"),signal.cex=1, 
       file="jpg",file.name="20250116_metabolism_for_legend",dpi=300,file.output=TRUE,verbose=TRUE,
       points.alpha=225,legend.ncol=3, legend.pos="middle")

CMplot(df_plot, plot.type="m",col=c("gray50"),multraits=TRUE,threshold=3.378378e-10,threshold.lty=1,
       threshold.lwd=c(1,1), threshold.col=c("black","grey"),amplify=TRUE,
       chr.den.col=NULL, signal.col=c("#011627","#F71735","#41EAD4","#6883BA","#FF9F1C","#118AB2","#06D6A0"),signal.cex=1, 
       file="pdf",file.name="20250116_metabolism_for_legend",dpi=300,file.output=TRUE,verbose=TRUE,
       points.alpha=225,legend.ncol=3, legend.pos="middle")

CMplot(df_plot,plot.type="q",col=c("#011627","#F71735","#41EAD4","#6883BA","#FF9F1C","#118AB2","#06D6A0"),multraits=TRUE,ylab.pos=2,signal.pch=c(19,6,4),signal.cex=1.2,signal.col="red",
       conf.int=TRUE,box=FALSE,axis.cex=1,file="jpg",file.name="20250116_metabolism_qq_for_legend",dpi=300,file.output=TRUE,
       verbose=TRUE,width=8,height=8)

CMplot(df_plot,plot.type="q",col=c("#011627","#F71735","#41EAD4","#6883BA","#FF9F1C","#118AB2","#06D6A0"),multraits=TRUE,ylab.pos=2,signal.pch=c(19,6,4),signal.cex=1.2,signal.col="red",
       conf.int=TRUE,box=FALSE,axis.cex=1,file="pdf",file.name="20250116_metabolism_qq_for_legend",dpi=300,file.output=TRUE,
       verbose=TRUE,width=8,height=8)

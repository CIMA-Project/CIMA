#!bin/bash

gcta=/media/scPBMC1_AnalysisDisk1/huangzhuoli/Script_HPC/software_gaoyue/GCTA/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static
bfile=/CIMA/genetics/qc/for.gcta.mlma_2025/378samples
grm=/CIMA/genetics/qc/for.gcta.mlma_2025/378samplesgrm
qcovar=/CIMA/Data/GWAS/qcovar_meta_lipid.txt
covar=/CIMA/Data/GWAS/covar_meta.txt

METABOLIC=$1
pheno=/CIMA/Data/GWAS/meta/META_${METABOLIC}.txt
out=/CIMA/CIMA_r1/Result/GWAS_meta/META_${METABOLIC}_378

${gcta} --mlma --bfile ${bfile} --grm ${grm} --qcovar ${qcovar} --covar ${covar} --pheno ${pheno} --thread-num 1 --out ${out} && \
echo "** META_${METABOLIC} DONE **"

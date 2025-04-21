#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
os.environ["R_HOME"] = "/home/huangzhuoli/mambaforge/envs/tensorqtl/lib/R"
os.environ["R_LIBS_USER"] = "/home/huangzhuoli/mambaforge/envs/tensorqtl/lib/R/library"
import pandas as pd
import torch
import tensorqtl
from tensorqtl import  trans,rfunc
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"torch: {torch.__version__} (CUDA {torch.version.cuda}), device: {device}")
print(f"pandas: {pd.__version__}")


# In[2]:


import argparse
# 定义仅需要 'cell' 参数的解析器
parser = argparse.ArgumentParser(description='Define arguments')
parser.add_argument('--cell', type=str, required=True)

# 解析命令行参数
args = parser.parse_args()

# 获取并打印 'cell' 参数
cell = args.cell
print(f'processing_{cell}')


# In[6]:


#read_genotype_df and vatiant_df
genotype_df = pd.read_parquet('/CIMA/Data/413_sample_genotype.parquet')
variant_df = pd.read_parquet('/CIMA/Data/413_sample_variant_df.parquet')


# In[7]:


#筛选需要检测的位点
lead_esnp_all_df =  pd.read_csv('/CIMA/Result/20250108_cis_caQTL_studywise_sig.csv',index_col=0) 

lead_esnp_all_df = lead_esnp_all_df[lead_esnp_all_df["celltype"] == cell]

lead_esnp_celltype = lead_esnp_all_df['variant_id'].unique()

genotype_df = genotype_df.loc[genotype_df.index.isin(lead_esnp_celltype),]

variant_df = variant_df.loc[variant_df.index.isin(lead_esnp_celltype),]


# In[8]:


#read_expression_file and cov
expression_bed = f'/CIMA/Data/caQTL/bed_file/{cell}.bed'
covariates_file = f'/CIMA/Data/caQTL/peer/peer_rez/factor/{cell}.csv'


# In[9]:


# load phenotypes and covariates
# align phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep=',', index_col=0)
phenotype_df = phenotype_df[covariates_df.index]


# In[10]:


# decited the number of peer factor
sample_size = covariates_df.shape[0]
if sample_size <= 150:
    pf = list('pf' + str(i) for i in range(1, 16))
elif sample_size <= 250:
    pf = list('pf' + str(i) for i in range(1, 31))
elif sample_size > 250:
    pf = list('pf' + str(i) for i in range(1, 36))

#pf = list('pf' + str(i) for i in range(1, 36))
info = list(['Age', 'sex', 'PC1', 'PC2'])
cova_need = info + pf
covariates_df = covariates_df.loc[:,cova_need]  


# In[11]:


trans_df = trans.map_trans(genotype_df, phenotype_df, covariates_df, batch_size=5000,
                           return_sparse=True, pval_threshold=1, maf_threshold=0.1)


# In[13]:
trans_df['variant_chr'] = trans_df['variant_id'].str.split('_').str[0].str.replace('chr', '')

trans_df['phenotype_chr'] = phenotype_pos_df.loc[trans_df['phenotype_id'],'chr'].values

trans_df['celltype'] = cell

trans_df.to_parquet(f'/CIMA/Result/caQTL_trans/{cell}.parquet')


# In[14]:

print("success")
print(f'{cell}_finish')


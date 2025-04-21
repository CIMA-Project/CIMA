#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
os.environ["R_HOME"] = "/home/huangzhuoli/mambaforge/envs/tensorqtl/lib/R"
os.environ["R_LIBS_USER"] = "/home/huangzhuoli/mambaforge/envs/tensorqtl/lib/R/library"
import pandas as pd
import torch
import tensorqtl
from tensorqtl import  cis
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"torch: {torch.__version__} (CUDA {torch.version.cuda}), device: {device}")
print(f"pandas: {pd.__version__}")


# In[3]:


import argparse
# 定义仅需要 'cell' 参数的解析器
parser = argparse.ArgumentParser(description='Define arguments')
parser.add_argument('--cell', type=str, required=True)

# 解析命令行参数
args = parser.parse_args()

# 获取并打印 'cell' 参数
cell = args.cell
print(f'processing_{cell}')


# In[ ]:


## Use tensorqtl
nominal_out_dir = f'/CIMA/Result/caQTL_tensorqtl/{cell}/'
if not os.path.exists(nominal_out_dir):
    os.makedirs(nominal_out_dir)
    print('dir created')


# In[ ]:


#read_genotype_df and vatiant_df
genotype_df = pd.read_parquet('/CIMA/Data/413_sample_genotype.parquet')
variant_df = pd.read_parquet('/CIMA/Data/413_sample_variant_df.parquet')


# In[ ]:


#read_expression_file and cov
expression_bed = f'/CIMA/Data/caQTL/bed_file/{cell}.bed'
covariates_file = f'/CIMA/Data/caQTL/peer/peer_rez/factor/{cell}.csv'


# In[ ]:


# load phenotypes and covariates
# align phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep=',', index_col=0)
phenotype_df = phenotype_df[covariates_df.index]


# In[ ]:


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


# In[ ]:


#run_cis_permutation
cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df, maf_threshold=0.1, seed=8888)


# In[ ]:

#后续统一再在jupyter notebook里面算
#tensorqtl.calculate_qvalues(cis_df)


# In[ ]:


cis_df['celltype'] = cell


# In[ ]:


cis_df.to_csv(f'{nominal_out_dir}all_lead_perm.csv')


# In[ ]:


cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, 'all', covariates_df, maf_threshold=0.1, output_dir=nominal_out_dir)

print(f'{cell}_done')

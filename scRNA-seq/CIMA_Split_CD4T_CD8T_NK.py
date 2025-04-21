import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('/media/AnalysisFastDisk/NatualCohort_T&NK_Cluster1st_5.2M.h5ad')

CD4T_temp1 = pd.read_csv('./CD4T_temp1_metadata.csv',index=0)
CD8T_temp1 = pd.read_csv('./CD8T_temp1_metadata.csv',index=0)
CD8T_temp2 = pd.read_csv('./CD8T_temp2_metadata.csv',index=0)
CD8T_temp3 = pd.read_csv('./CD8T_temp3_metadata.csv',index=0)

metadata = [CD4T_temp1, CD8T_temp1, CD8T_temp2, CD8T_temp3]
metadata = pd.concat(metadata)

metadata['celltype_temp'] = metadata['celltype_temp4']
metadata['celltype_temp'].fillna(metadata['celltype_temp1'],inplace=True)
metadata['celltype_temp'].fillna(metadata['celltype_temp2'],inplace=True)
metadata['celltype_temp'].fillna(metadata['celltype_temp3'],inplace=True)

adata.obs['celltype_TNK_1st'] = adata.obs['leiden_r2_n2']
adata.obs['celltype_TNK_1st'] = adata.obs.index.map(metadata['celltype_temp'].to_dict())
adata.obs['celltype_TNK_1st'].fillna(adata.obs['leiden_r2_n2'],inplace=True)

adata_NK = adata[adata.obs['celltype_TNK_1st'].isin(['NK','Cycling_NK','12','19','20'])].copy()
adata_NK.write('/media/AnalysisFastDisk/NatualCohort_NK_0.5M.h5ad')

adata_CD4T = adata[adata.obs['celltype_TNK_1st'].isin(['CD4T'])].copy()
adata_CD4T.write('/media/AnalysisFastDisk/NatualCohort_CD4T_2.3M.h5ad')

adata_CD8T = adata[~(adata.obs['celltype_TNK_1st'].isin(['NK','Cycling_NK','12','19','20','CD4T']))].copy()
adata_CD8T.write('/media/AnalysisFastDisk/NatualCohort_CD8T_2.3M.h5ad')
import anndata as ad
import pandas as pd
import numpy as np
from tqdm import tqdm

# Import Metacell Data
rna = ad.read_h5ad("CIMA/scRNA/CIMA_scRNA_Metacell.h5ad")
atac = ad.read_h5ad("CIMA/scATAC/CIMA_scATAC_PeakRecalling_Metacell.h5ad")

rna = rna[rna.obs['cell_type_l4'].isin(list(atac.obs['cell_type_l4'].unique()))]

new_barcodes = []
rna.obs['new_barcode'] = ''
atac.obs['new_barcode'] = ''

# Generating Pseudo-multi-omics Metacell Data
for (sample, celltype), rna_group in tqdm(rna.obs.groupby(['sample', 'cell_type_l4'])):
    atac_group = atac.obs[(atac.obs['sample'] == sample) & (atac.obs['cell_type_l4'] == celltype)]
    
    if not rna_group.empty and not atac_group.empty:
        
        min_cell_count = min(len(rna_group), len(atac_group))
        
        selected_rna_cells = rna_group.sample(n=min_cell_count, random_state=42)
        selected_atac_cells = atac_group.sample(n=min_cell_count, random_state=42)
        
        for i in range(min_cell_count):
            new_barcode = f"{sample}_{celltype}_{i+1}"
            
            rna.obs.at[selected_rna_cells.index[i], 'new_barcode'] = new_barcode
            rna.obs.at[selected_atac_cells.index[i], 'new_barcode'] = new_barcode
            
            new_barcodes.append({
                'new_barcode': new_barcode,
                'sample': sample,
                'cell_type': celltype
            })

barcode_info_df = pd.DataFrame(new_barcodes)
barcode_info_df.to_csv("CIMA/Pseudo_multiomics/Pseudo_multiomics_barcode_info.csv", index=False)

rna_filtered = rna[rna.obs['new_barcode'] != ''].copy()
rna_filtered.write_h5ad('CIMA/Pseudo_multiomics/CIMA_scRNA_Pseudo_multiomics.h5ad')

atac_filtered = atac[atac.obs['new_barcode'] != ''].copy()
atac_filtered.write_h5ad('CIMA/Pseudo_multiomics/CIMA_scATAC_Pseudo_multiomics.h5ad')
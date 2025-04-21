import warnings
import os
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import dill
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

_stderr = sys.stderr
null = open(os.devnull,'wb')
work_dir = 'CIMA/scenicplus/'

scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus_pseudo_multiomics_obj.pkl'), 'rb'))

from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons
apply_std_filtering_to_eRegulons(scplus_obj)

scplus_obj.uns['eRegulon_metadata_filtered'].to_csv(os.path.join(work_dir, 'All_eRegulon_metadata_filtered.csv'))


#eRegulon enrichment scores
from scenicplus.eregulon_enrichment import score_eRegulons
region_ranking = dill.load(open(os.path.join(work_dir, 'scenicplus_custom_cell_type_l4/region_ranking.pkl'), 'rb')) #load ranking calculated using the wrapper function
gene_ranking = dill.load(open(os.path.join(work_dir, 'scenicplus_custom_cell_type_l4/gene_ranking.pkl'), 'rb')) #load ranking calculated using the wrapper function

score_eRegulons(scplus_obj,
                ranking = region_ranking,
                eRegulon_signatures_key = 'eRegulon_signatures_filtered',
                key_added = 'eRegulon_AUC_filtered',
                enrichment_type= 'region',
                auc_threshold = 0.05,
                normalize = False,
                n_cpu = 1)

score_eRegulons(scplus_obj,
                gene_ranking,
                eRegulon_signatures_key = 'eRegulon_signatures_filtered',
                key_added = 'eRegulon_AUC_filtered',
                enrichment_type = 'gene',
                auc_threshold = 0.05,
                normalize= False,
                n_cpu = 1)

#eRegulon dimensionality reduction
from scenicplus.dimensionality_reduction import run_eRegulons_tsne, run_eRegulons_umap
run_eRegulons_tsne(
    scplus_obj = scplus_obj,
    auc_key = 'eRegulon_AUC_filtered',
    reduction_name = 'eRegulons_tSNE',random_state=555 #overwrite previously calculated tSNE
)
run_eRegulons_umap(
    scplus_obj = scplus_obj,
    auc_key = 'eRegulon_AUC_filtered',
    reduction_name = 'eRegulons_UMAP',min_dist=0.2,random_state=123 #overwrite previously calculated UMAP
)


from scenicplus.cistromes import TF_cistrome_correlation, generate_pseudobulks

generate_pseudobulks(
        scplus_obj = scplus_obj,
        nr_cells = 10,
        variable = 'GEX_cell_type_l4',
        auc_key = 'eRegulon_AUC_filtered',
        signature_key = 'Gene_based')


# TF cistrome correlation in each celltype l1

cell_types = {
    'B': 'B',
    'Myeloid': 'Myeloid',
    'CD4T': 'CD4',
    'CD8T': 'CD8',
    'NK': 'NK'
}


thresholds = {
    'rho': [-0.8, 0.8],
    'n_targets': 0
}


for cell_type, acc_key in cell_types.items():

    celltype_l4 = list(scplus_obj.metadata_cell[scplus_obj.metadata_cell['ACC_cell_type_l1'] == acc_key]['GEX_cell_type_l4'].unique())
    

    TF_cistrome_correlation(
        scplus_obj,
        use_pseudobulk=True,
        variable='GEX_cell_type_l4',
        auc_key='eRegulon_AUC_filtered',
        signature_key='Gene_based',
        out_key=f'filtered_gene_based_{cell_type}',
        subset=celltype_l4
    )
    

    pval_key = f'filtered_gene_based_{cell_type}'
    scplus_obj.uns['TF_cistrome_correlation'][pval_key]['Adjusted_p-value_nonzero'] = (
        scplus_obj.uns['TF_cistrome_correlation'][pval_key]['Adjusted_p-value'].replace({0.000000e+00: 1e-300}, inplace=False)
    )
    

    scplus_obj.uns['TF_cistrome_correlation'][pval_key].to_csv(
        os.path.join(work_dir, f'TF_cistrome_correlation_filtered_gene_based_{cell_type}.csv')
    )
    

    n_targets = [int(x.split('(')[1].replace('g)', '')) for x in scplus_obj.uns['TF_cistrome_correlation'][pval_key]['Cistrome']]
    rho = scplus_obj.uns['TF_cistrome_correlation'][pval_key]['Rho'].to_list()
    adj_pval = scplus_obj.uns['TF_cistrome_correlation'][pval_key]['Adjusted_p-value_nonzero'].to_list()
    

    fig, ax = plt.subplots(figsize=(10, 5))
    sc = ax.scatter(rho, n_targets, c=-np.log10(adj_pval), s=10)
    ax.set_xlabel('Correlation coefficient')
    ax.set_ylabel('nr. target regions')
    ax.vlines(x=thresholds['rho'], ymin=0, ymax=max(n_targets), color='black', ls='dashed', lw=1)
    ax.text(x=thresholds['rho'][0], y=max(n_targets), s=str(thresholds['rho'][0]))
    ax.text(x=thresholds['rho'][1], y=max(n_targets), s=str(thresholds['rho'][1]))
    sns.despine(ax=ax)
    fig.colorbar(sc, label='-log10(adjusted_pvalue)', ax=ax)
    plt.show()
    

    selected_cistromes = scplus_obj.uns['TF_cistrome_correlation'][pval_key].loc[
        np.logical_or(
            scplus_obj.uns['TF_cistrome_correlation'][pval_key]['Rho'] > thresholds['rho'][1],
            scplus_obj.uns['TF_cistrome_correlation'][pval_key]['Rho'] < thresholds['rho'][0]
        )
    ]['Cistrome'].to_list()
    

    selected_eRegulons = [x.split('_(')[0] for x in selected_cistromes]
    selected_eRegulons_gene_sig = [
        x for x in scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'].keys()
        if x.split('_(')[0] in selected_eRegulons
    ]
    selected_eRegulons_region_sig = [
        x for x in scplus_obj.uns['eRegulon_signatures_filtered']['Region_based'].keys()
        if x.split('_(')[0] in selected_eRegulons
    ]

    scplus_obj.uns[f'selected_{cell_type}_eRegulon'] = {
        'Gene_based': selected_eRegulons_gene_sig,
        'Region_based': selected_eRegulons_region_sig
    }
    print(f'Selected {len(selected_eRegulons_gene_sig)} eRegulons for {cell_type}')
    

    with open(os.path.join(work_dir, f'selected_eRegulons_{cell_type}_0.8.txt'), 'w') as f:
        for eRegulon_name in scplus_obj.uns[f'selected_{cell_type}_eRegulon']['Gene_based']:
            f.write(eRegulon_name.rsplit('_', 1)[0] + '\n')

# TF cistrome correlation in all celltypes
TF_cistrome_correlation(
            scplus_obj,
            use_pseudobulk = True,
            variable = 'GEX_cell_type_l4',
            auc_key = 'eRegulon_AUC_filtered',
            signature_key = 'Gene_based',
            out_key = 'filtered_gene_based_All')


pval_key = f'filtered_gene_based_All'
scplus_obj.uns['TF_cistrome_correlation'][pval_key]['Adjusted_p-value_nonzero'] = (
    scplus_obj.uns['TF_cistrome_correlation'][pval_key]['Adjusted_p-value'].replace({0.000000e+00: 1e-300}, inplace=False)
)


scplus_obj.uns['TF_cistrome_correlation'][pval_key].to_csv(
    os.path.join(work_dir, f'TF_cistrome_correlation_filtered_gene_based_All.csv')
)


n_targets = [int(x.split('(')[1].replace('g)', '')) for x in scplus_obj.uns['TF_cistrome_correlation'][pval_key]['Cistrome']]
rho = scplus_obj.uns['TF_cistrome_correlation'][pval_key]['Rho'].to_list()
adj_pval = scplus_obj.uns['TF_cistrome_correlation'][pval_key]['Adjusted_p-value_nonzero'].to_list()


thresholds = {
    'rho': [-0.6, 0.6],
    'n_targets': 0
}


fig, ax = plt.subplots(figsize=(10, 5))
sc = ax.scatter(rho, n_targets, c=-np.log10(adj_pval), s=10)
ax.set_xlabel('Correlation coefficient')
ax.set_ylabel('nr. target regions')
ax.vlines(x=thresholds['rho'], ymin=0, ymax=max(n_targets), color='black', ls='dashed', lw=1)
ax.text(x=thresholds['rho'][0], y=max(n_targets), s=str(thresholds['rho'][0]))
ax.text(x=thresholds['rho'][1], y=max(n_targets), s=str(thresholds['rho'][1]))
sns.despine(ax=ax)
fig.colorbar(sc, label='-log10(adjusted_pvalue)', ax=ax)
plt.show()


selected_cistromes = scplus_obj.uns['TF_cistrome_correlation'][pval_key].loc[
    np.logical_or(
        scplus_obj.uns['TF_cistrome_correlation'][pval_key]['Rho'] > thresholds['rho'][1],
        scplus_obj.uns['TF_cistrome_correlation'][pval_key]['Rho'] < thresholds['rho'][0]
    )
]['Cistrome'].to_list()


selected_eRegulons = [x.split('_(')[0] for x in selected_cistromes]
selected_eRegulons_gene_sig = [
    x for x in scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'].keys()
    if x.split('_(')[0] in selected_eRegulons
]
selected_eRegulons_region_sig = [
    x for x in scplus_obj.uns['eRegulon_signatures_filtered']['Region_based'].keys()
    if x.split('_(')[0] in selected_eRegulons
]


scplus_obj.uns[f'selected_All_eRegulon'] = {
    'Gene_based': selected_eRegulons_gene_sig,
    'Region_based': selected_eRegulons_region_sig
}
print(f'Selected {len(selected_eRegulons_gene_sig)} eRegulons for All')


with open(os.path.join(work_dir, f'selected_eRegulons_All_0.6.txt'), 'w') as f:
    for eRegulon_name in scplus_obj.uns[f'selected_All_eRegulon']['Gene_based']:
        f.write(eRegulon_name.rsplit('_', 1)[0] + '\n')


# Export loom file
from scenicplus.eregulon_enrichment import binarize_AUC
binarize_AUC(scplus_obj,
             auc_key='eRegulon_AUC_filtered',
             out_key='eRegulon_AUC_filtered_thresholds',
             signature_keys=['Gene_based'],
             n_cpu=96)

from scenicplus.loom import *
export_to_loom(scplus_obj,
       signature_key = 'Gene_based',
       eRegulon_metadata_key = 'eRegulon_metadata_filtered',
       auc_key = 'eRegulon_AUC_filtered',
       auc_thr_key = 'eRegulon_AUC_filtered_thresholds',
       keep_direct_and_extended_if_not_direct = True,
       tree_structure = (),
       title =  'Gene based eGRN',
       nomenclature = 'hg38',
       out_fname=os.path.join(work_dir,'SCENIC+_gene_based_filter.loom'))

# Peak to gene file
from scenicplus.enhancer_to_gene import export_to_UCSC_interact
r2g_data = export_to_UCSC_interact(scplus_obj,
                            'hsapiens',
                            os.path.join(work_dir,'CIMA_cell_type_l4_r2g.rho.bed'),
                            path_bedToBigBed='./',
                            bigbed_outfile=os.path.join(work_dir,'CIMA_cell_type_l4_r2g.rho.interact.bb'),
                            region_to_gene_key='region_to_gene',
                            pbm_host="http://sep2019.archive.ensembl.org/",
                            assembly='hg38',
                            ucsc_track_name='cell_type_l4_R2G',
                            ucsc_description='Region_to_gene in the CIMA celltype l4',
                            cmap_neg='Reds',
                            cmap_pos='Greens',
                            key_for_color='rho',
                            vmin=-1,
                            vmax=1,
                            scale_by_gene=False,
                            subset_for_eRegulons_regions=True,
                            eRegulons_key='eRegulons')

# AUC correlation heatmap

combined_selected_eRegulon_list = (
    scplus_obj.uns['selected_All_eRegulon']['Gene_based'] +
    scplus_obj.uns['selected_B_eRegulon']['Gene_based'] +
    scplus_obj.uns['selected_CD4T_eRegulon']['Gene_based'] +
    scplus_obj.uns['selected_CD8T_eRegulon']['Gene_based'] +
    scplus_obj.uns['selected_NK_eRegulon']['Gene_based'] +
    scplus_obj.uns['selected_Myeloid_eRegulon']['Gene_based']
)

filtered_list = [tf for tf in combined_selected_eRegulon_list if '_+_' in tf]
unique_filtered_list = list(set(filtered_list))
from scenicplus.plotting.correlation_plot import *
correlation_heatmap(scplus_obj,
                    auc_key = 'eRegulon_AUC_filtered',
                    signature_keys = ['Gene_based'],
                    scale=False,
                    selected_regulons = unique_filtered_list,
                    fcluster_threshold = 0.1,use_plotly='seaborn',cmap='GnBu',
                    fontsize = 5,
                    save = os.path.join(work_dir,'eRegulon_AUC_filtered_Gene_based_correlation_heatmap.pdf')
                   )
correlation_heatmap(scplus_obj,
                    auc_key = 'eRegulon_AUC_filtered',
                    signature_keys = ['Region_based'],
                    scale=False,
                    selected_regulons = unique_filtered_list,
                    fcluster_threshold = 0.1,use_plotly='seaborn',cmap='GnBu',
                    fontsize = 5,
                    save = os.path.join(work_dir,'eRegulon_AUC_filtered_Region_based_correlation_heatmap.pdf')
                   )

# Export AUC Matrix
auc_df = scplus_obj.uns['eRegulon_AUC_filtered']['Gene_based']
auc_df.to_parquet(os.path.join(work_dir,'scenicplus_celltype_l4_AUC.parquet'))

scplus_obj.metadata_cell.to_parquet(os.path.join(work_dir,'scenicplus_celltype_l4_metadata.parquet'))

auc_df_with_metadata = pd.merge(
    auc_df,
    scplus_obj.metadata_cell[['ACC_sample','ACC_cell_type_l4','ACC_cell_type_l3','ACC_cell_type_l2','ACC_cell_type_l1','ACC_sex','ACC_age','ACC_age_group','ACC_cell_type_l4_age_sex_group']],
    left_index=True,
    right_index=True)
auc_df_with_metadata.to_parquet(os.path.join(work_dir,'scenicplus_celltype_l4_metadata_with_Metadata.parquet'))
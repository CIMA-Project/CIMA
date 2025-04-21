import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

auc_df = pd.read_parquet('CIMA/scenicplus/scenicplus_celltype_l4_AUC.parquet')
adata = sc.read_h5ad('CIMA/Pseudo_multiomics/CIMA_scRNA_Pseudo_multiomics.h5ad')
adata.obs.index = adata.obs['new_barcode']

Mono_adata = adata[adata.obs['cell_type_l3'].isin(['cMono','ncMono','intMono'])]
Mono_adata.layers['counts'] = Mono_adata.X.copy()
sc.pp.normalize_total(Mono_adata, target_sum=1e4)
sc.pp.log1p(Mono_adata)
sc.pp.highly_variable_genes(Mono_adata, flavor='seurat_v3', n_top_genes=1000)
Mono_adata = Mono_adata[:, Mono_adata.var.highly_variable]
sc.pp.scale(Mono_adata, max_value=10)
sc.tl.pca(Mono_adata, svd_solver="arpack")
sc.pl.pca_variance_ratio(Mono_adata, log=True)
sc.external.pp.harmony_integrate(Mono_adata,key='sample')
sc.pp.neighbors(Mono_adata, n_neighbors=30, n_pcs=5,use_rep='X_pca_harmony')
sc.tl.umap(Mono_adata)

Mono_adata.X = Mono_adata.layers['counts'].copy()
Mono_adata.X = Mono_adata.X.astype('float32')
sc.pp.calculate_qc_metrics(Mono_adata, inplace=True)

tnode = sct.train.Trainer(Mono_adata, loss_mode='nb', alpha_recon_lec=0.5, alpha_recon_lode=0.5)
tnode.train()
Mono_adata.obs['ptime'] = tnode.get_time()

mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)
Mono_adata.obsm['X_TNODE'] = mix_zs
Mono_adata.obsm['X_VF'] = tnode.get_vector_field(Mono_adata.obs['ptime'].values, Mono_adata.obsm['X_TNODE'])

Mono_adata.obs = pd.merge(Mono_adata.obs,auc_df,left_index=True,right_index=True)

sct.vf.plot_vector_field(Mono_adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE',
                         color='cell_type_l4', legend_loc='none', frameon=False, size=100, alpha=0.5,reverse=True,
                         palette=cell_type_l4_colors, save='CIMA/scTour/scTour_Mono_cell_type_l4_TF_UMAP/scTour_Mono_cell_type_l4_Umap.pdf')

for TF in list(Mono_adata.obs.columns[13:]):
    sct.vf.plot_vector_field(Mono_adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE',
                             color=TF, legend_loc='none', cmap = plt.cm.OrRd, frameon=False, size=100, alpha=1, reverse=True,
                             save='CIMA/scTour/scTour_Mono_cell_type_l4_TF_UMAP/scTour_Mono_cell_type_l4_TF_'+TF+'.pdf')

# Pseudotime Distribution by cell type l4
median_pseudotime = Mono_adata.obs.groupby('cell_type_l4')['ptime'].median()
sorted_cell_types = median_pseudotime.sort_values().index
plt.figure(figsize=(12, 8))
sns.boxplot(x='cell_type_l4', y='ptime', data=Mono_adata.obs, order=sorted_cell_types)
plt.title('Pseudotime Distribution by Cell Type')
plt.xlabel('Cell Type l4')
plt.ylabel('Pseudotime')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()


sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
df = Mono_adata.obs[['cell_type_l4', 'ptime']]
median_pseudotime = df.groupby('cell_type_l4')['ptime'].median().sort_values()
df['cell_type_l4'] = pd.Categorical(df['cell_type_l4'], categories=median_pseudotime.index, ordered=True)
pal = sns.cubehelix_palette(len(median_pseudotime), rot=-.25, light=.7)

g = sns.FacetGrid(df, row="cell_type_l4", hue="cell_type_l4", aspect=15, height=.8, palette=pal)
g.map(sns.kdeplot, "ptime",
      bw_adjust=.5, clip_on=False,
      fill=True, alpha=1, linewidth=1)
g.map(sns.kdeplot, "ptime", clip_on=False, color="w", lw=2, bw_adjust=.5)
g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

def label(x, color, label):
    ax = plt.gca()
    ax.text(0, .2, label, fontweight="bold", color=color,
            ha="left", va="center", transform=ax.transAxes)
g.map(label, "ptime")

g.figure.subplots_adjust(hspace=-.25)
g.set_titles("")
g.set(yticks=[], ylabel="")
g.despine(bottom=True, left=True)
g.savefig("Mono_Pseudotime_Distribution_by_cell_type_l4.pdf")
plt.show()

Mono_adata.write('CIMA/scTour/Mono_Pseudo_Multiomics_Pseudotim.h5ad')


# Calculate the correlation between transcription factor activity in monocytes and pseudotime
rss_df = pd.read_csv('rel_data_Activators_Top15.csv',index_col=0)

cell_type_list = adata.obs['cell_type_l4'].unique()
filtered_df = rss_df[(rss_df['Top_RSS'] == True) & (rss_df['Cell_type'].isin(cell_type_list))]
filtered_tf_list = filtered_df['Regulon'].unique()

ptime = Mono_adata.obs['ptime']
tf_activity = Mono_adata.obs.iloc[:, 13:]
ptime_reversed = 1 - ptime

correlations = []
p_values = []
for tf in filtered_tf_list:
    corr, p_val = pearsonr(tf_activity[tf], ptime_reversed)
    correlations.append(corr)
    p_values.append(p_val)

correlation_df = pd.DataFrame({
    'Transcription Factor': filtered_tf_list,
    'Correlation': correlations,
    'P-value': p_values
})

significant_df = correlation_df[correlation_df['P-value'] < 0.05]
significant_df = significant_df.sort_values(by='Correlation', ascending=False)
colors = plt.cm.coolwarm((significant_df['Correlation'] - significant_df['Correlation'].min()) / 
                         (significant_df['Correlation'].max() - significant_df['Correlation'].min()))

fig, ax = plt.subplots(figsize=(12, 6))
for i, (tf, row) in enumerate(significant_df.iterrows()):
    color = colors[i]
    ax.plot(i, row['Correlation'], marker='o', linestyle='-', color=color, markersize=10)
ax.set_xticks(range(len(significant_df)))
ax.set_xticklabels(significant_df['Transcription Factor'], rotation=90)
ax.set_xlabel('Transcription Factors')
ax.set_ylabel('Correlation with Pseudotime (Reversed)')
ax.set_title('Significant Correlations between Pseudotime and Transcription Factors')

sm = plt.cm.ScalarMappable(cmap='coolwarm', 
                           norm=plt.Normalize(vmin=significant_df['Correlation'].min(), 
                                              vmax=significant_df['Correlation'].max()))
sm.set_array([])

cbar = fig.colorbar(sm, ax=ax)
cbar.set_label('Correlation')

plt.tight_layout()
plt.show()
# Calculation of correlation between TF activity and Age
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

def analyze_celltype_correlation(
    AUC_df,
    RSS_df,
    parent_celltype,
    output_prefix,
    parent_celltype_col='ACC_cell_type_l1',
    cell_type_col='ACC_cell_type_l4',
    value_col='ACC_age',
    metadata_columns=8,
    alpha=0.05
):
    """
    Analyzing transcription factor-age correlations for each daughter cell type under a given parent cell type
    
    Parameters:
    AUC_df -- DataFrame containing AUC values and metadata
    RSS_df -- DataFrame containing RSS filter results
    subset_celltype -- Subset cell type list to be analyzed (e.g. ['CD8'] or ['cMono','ncMono'])
    output_prefix -- prefix of the output file
    subset_celltype_col -- name of the column where the subset cell type is located
    cell_type_col -- the column name of the cell type level used to calculate the correlation
    value_col -- Column name of the value used to calculate the correlation
    metadata_columns -- number of metadata columns (last n columns)
    alpha -- significance threshold
    """
    
    AUC_sub_df = AUC_df[AUC_df[parent_celltype_col].isin(parent_celltype)]
    cell_type_list = AUC_sub_df[cell_type_col].unique().tolist()
    
    # Filter RSS Data
    filtered_df = RSS_df[(RSS_df['Top_RSS'] == True) & 
                        (RSS_df['Cell_type'].isin(cell_type_list))]
    filtered_tf_list = filtered_df['Regulon'].unique()
    
    # Get TF columns (exclude last metadata column)
    transcription_factors = AUC_sub_df.columns[:-metadata_columns]
    transcription_factors_filtered = transcription_factors[
        transcription_factors.str.startswith(tuple(filtered_tf_list))
    ]
    
    correlation_df = pd.DataFrame(index=transcription_factors_filtered,
                                 columns=AUC_sub_df[cell_type_col].unique())
    p_value_df = pd.DataFrame(index=transcription_factors_filtered,
                            columns=AUC_sub_df[cell_type_col].unique())
    
    # Calculate correlation
    for cell_type in AUC_sub_df[cell_type_col].unique():
        subset = AUC_sub_df[AUC_sub_df[cell_type_col] == cell_type]
        for tf in transcription_factors_filtered:
            if len(subset) >= 2:
                corr, p_value = pearsonr(subset[tf], subset[value_col])
            else:
                corr = p_value = np.nan
            correlation_df.loc[tf, cell_type] = corr
            p_value_df.loc[tf, cell_type] = p_value
    
    correlation_df = correlation_df.astype(float)
    correlation_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    correlation_df.dropna(how='all', inplace=True)
    
    def plot_clustermap(corr_df, pval_df, filename):
        g = sns.clustermap(corr_df, cmap="Spectral_r", vmin=-0.4, vmax=0.4,
                          linewidths=0, figsize=(8,17))
        ax = g.ax_heatmap

        for i, tf in enumerate(corr_df.index):
            for j, ct in enumerate(corr_df.columns):
                if pval_df.loc[tf, ct] < alpha:
                    ax.text(j + 0.5, i + 0.5, '*', 
                           ha='center', va='center',
                           color='black', fontsize=10)
        
        plt.title(f'Correlation Analysis: {parent_celltype} ~ {value_col}')
        plt.savefig(f'{output_prefix}_Correlation_TF_{value_col}.pdf', bbox_inches='tight')
        plt.close()
    
    plot_clustermap(correlation_df, p_value_df, output_prefix)
    
    return correlation_df, p_value_df


AUC_df = pd.read_parquet('scenicplus_celltype_l4_metadata_with_Metadata.parquet')
RSS_df = pd.read_csv('rel_data_Activators_Top15.csv', index_col=0)

# Figure 3I
CD8_corr, CD8_pval = analyze_celltype_correlation(
    AUC_df, RSS_df, ['CD8'], 'CD8T',
    'ACC_cell_type_l1', 'ACC_cell_type_l4','ACC_age'
    )
CD8_corr.to_csv('CD8T_RSS_top15_TF_correlation_with_Age.csv')
CD8_pval.to_csv('CD8T_RSS_top15_TF_correlation_with_Age_Pvalue.csv')


# Figure 3J
Mono_corr, Mono_pval = analyze_celltype_correlation(
    AUC_df, RSS_df, ['cMono','intMono','ncMono'], 'Mono',
    'ACC_cell_type_l3', 'ACC_cell_type_l4','ACC_age'
)
Mono_corr.to_csv('Mono_RSS_top15_TF_correlation_with_Age.csv')
Mono_pval.to_csv('Mono_RSS_top15_TF_correlation_with_Age_Pvalue.csv')

# Calculation of the correlation between transcription factor AUC values and age in all cell types
transcription_factors = AUC_df.columns[:-8] 

correlation_df = pd.DataFrame(index=transcription_factors, columns=AUC_df['ACC_cell_type_l4'].unique())
p_value_df = pd.DataFrame(index=transcription_factors, columns=AUC_df['ACC_cell_type_l4'].unique())

for cell_type in AUC_df['ACC_cell_type_l4'].unique():
    subset = AUC_df[AUC_df['ACC_cell_type_l4'] == cell_type]
    for tf in transcription_factors:
        corr, p_value = pearsonr(subset[tf], subset['ACC_age'])
        correlation_df.loc[tf, cell_type] = corr
        p_value_df.loc[tf, cell_type] = p_value

correlation_df = correlation_df.astype(float)

correlation_df.replace([np.inf, -np.inf], np.nan, inplace=True)
correlation_df.dropna(inplace=True)

correlation_df.to_csv('Cell_type_l4_TF_correlation_with_Age.csv')
p_value_df.to_csv('Cell_type_l4_TF_correlation_with_Age_Pvalue.csv')



# cell type-specific aging GRNs plot
import scanpy as sc
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'


def create_eregulon_network(
    cell_type, 
    tf_set,
    rna_data,
    peak_df,
    eRegulon_metadata,
    cell_type_col='cell_type_l4',
    tf_col='TF',
    region_col='Region',
    gene_col='Gene',
    group_col='Age_group',
    group_value='Old',
    deg_pval_thresh=0.05,
    deg_logfc_thresh=0.5,
    tf_prefix='TF-',
    output_path=None,
    figsize=(10, 8),
    dpi=300
):
    """
    Creating and visualizing eRegulon regulatory networks
    
    parameters:
    cell_type -- target cell type
    tf_set -- set of transcription factors to include
    rna_data -- single cell RNA data (AnnData format)
    peak_df -- DataFrame containing peak information
    eRegulon_metadata -- eRegulon metadata
    cell_type_col -- Cell type column name
    tf_col -- tf information column
    region_col -- region information column
    gene_col -- gene information column
    group_col -- column name of age group
    group_value -- value of age group to select up-regulated differentially expressed genes
    deg_pval_thresh -- differential gene p-value threshold
    deg_logfc_thresh -- Differential gene logFC thresholds
    tf_prefix -- TF node prefix
    output_path -- output file path
    figsize -- graph size
    dpi -- image resolution
    """
    
    # Get the differentially expressed genes of the specified cell type in the group
    cell_type_data = rna_data[rna_data.obs[cell_type_col] == cell_type].copy()
    sc.tl.rank_genes_groups(cell_type_data, group_col, method="wilcoxon")
    DEG_df = sc.get.rank_genes_groups_df(cell_type_data,None)
    DEG_sig_df = DEG_df[(DEG_df['group'] == group_value) & 
                       (DEG_df['pvals_adj'] < deg_pval_thresh) & 
                       (DEG_df['logfoldchanges'] > deg_logfc_thresh)]
    DEG_sig_list = set(DEG_sig_df['names'])
    
    # Get the open peak set
    peak_set = set(peak_df[peak_df[cell_type] == True].index)
    
    # Filtering eRegulon metadata
    filtered_eRegulon = eRegulon_metadata[
        (eRegulon_metadata[tf_col].isin(tf_set)) &
        (eRegulon_metadata[gene_col].isin(DEG_sig_list)) &
        (eRegulon_metadata[region_col].isin(peak_set))
    ].copy()
    filtered_eRegulon[tf_col] = tf_prefix + filtered_eRegulon[tf_col]
    
    # Creating Network
    G = nx.MultiDiGraph()
    for _, row in filtered_eRegulon.iterrows():
        tf_node = row[tf_col]
        region_node = row[region_col]
        gene_node = row[gene_col]
        
        G.add_node(tf_node, type='TF', label=row[tf_col])
        G.add_node(region_node, type='Region')
        G.add_node(gene_node, type='Gene', label=row[gene_col])
        
        G.add_edge(tf_node, region_node, type='TF-Region', tf=row[tf_col])
        G.add_edge(region_node, gene_node, type='Region-Gene', tf=row[tf_col])
    
    # Node Settings
    node_shapes = {'TF': 'h', 'Region': 'D', 'Gene': 'o'}
    node_colors = {'Region': 'blue'}
    node_sizes = {'TF': 1800, 'Region': 60, 'Gene': 300}
    node_type = nx.get_node_attributes(G, 'type')
    node_labels = nx.get_node_attributes(G, 'label')
    
    # Color Map Creating
    tf_labels = list(set(nx.get_edge_attributes(G, 'tf').values()))
    color_map = plt.cm.get_cmap('tab20', len(tf_labels))
    tf_color_map = {tf: color_map(i) for i, tf in enumerate(tf_labels)}
    
    # Hierarchical layout
    def layered_layout(G, layers):
        pos = {}
        for i, (layer, nodes) in enumerate(layers.items()):
            angle_step = 2 * np.pi / len(nodes)
            radius = (i + 1) * 1.5
            for j, node in enumerate(nodes):
                angle = j * angle_step
                pos[node] = (radius * np.cos(angle), radius * np.sin(angle))
        return pos
    
    # Node Layering
    node_type = nx.get_node_attributes(G, 'type')
    layers = {
        'TF': [n for n in G.nodes if node_type[n] == 'TF'],
        'Region': [n for n in G.nodes if node_type[n] == 'Region'],
        'Gene': [n for n in G.nodes if node_type[n] == 'Gene']
    }
    
    pos = layered_layout(G, layers)
    
    # Configure DEG node color
    vmin = DEG_sig_df['logfoldchanges'].min()
    vmax = DEG_sig_df['logfoldchanges'].max()
    norm = mcolors.TwoSlopeNorm(vcenter=0, vmax=vmax)
    cmap = plt.cm.bwr
    
    gene_colors = []
    for gene in layers['Gene']:
        gene_name = gene
        if gene_name in DEG_df['names'].values:
            logfc = DEG_df[DEG_df['names'] == gene_name]['logfoldchanges'].values[0]
            gene_colors.append(cmap(norm(logfc)))
        else:
            gene_colors.append('gray')
    
    plt.figure(figsize=figsize, dpi=dpi)
    
    nx.draw_networkx_nodes(
        G, pos, nodelist=layers['TF'],
        node_shape=node_shapes['TF'],
        node_size=node_sizes['TF'],
        node_color=[tf_color_map[node_labels.get(n, '')] for n in layers['TF']],
        label='TF'
    )
    
    nx.draw_networkx_nodes(
        G, pos, nodelist=layers['Region'],
        node_shape=node_shapes['Region'],
        node_size=node_sizes['Region'],
        node_color=node_colors['Region'],
        label='Region'
    )
    
    nx.draw_networkx_nodes(
        G, pos, nodelist=layers['Gene'],
        node_shape=node_shapes['Gene'],
        node_size=node_sizes['Gene'],
        node_color=gene_colors,
        label='Gene'
    )
    

    for i, (f, t, data) in enumerate(G.edges(data=True)):
        rad = 0.2 + 0.1*(i%3) if data['type'] == 'Region-Gene' else 0.2
        nx.draw_networkx_edges(
            G, pos, edgelist=[(f, t)],
            edge_color=[tf_color_map[data['tf']]],
            connectionstyle=f'arc3,rad={rad}'
        )
    

    node_labels = nx.get_node_attributes(G, 'label')
    labels_to_draw = {n: l for n, l in node_labels.items() if node_type[n] != 'Region'}
    nx.draw_networkx_labels(G, pos, labels=labels_to_draw, font_size=5)
    

    plt.title(f'{group_col}-{group_value} eRegulon Network in {cell_type}')
    plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), label='logfoldchanges')
    plt.legend(scatterpoints=1)
    

    if output_path:
        plt.savefig(output_path, bbox_inches='tight', format='pdf')
    plt.show()
    
    return G

eRegulon_metadata = pd.read_csv('CIMA/scenicplus/All_eRegulon_metadata_filtered.csv',index_col=0)
rna = sc.read_h5ad('CIMA/Pseudo_multiomics/CIMA_scRNA_Pseudo_multiomics.h5ad')
peak_df = pd.read_csv('CIMA/CIMA_Cell_type_l4_Peak.csv',index_col=0)

ncRNA_genes = list(rna.var_names[rna.var_names.str.match(r'^[A-Z][A-Z][0-9].*\.[0-9]')])
LINC_genes = list(rna.var_names[rna.var_names.str.match(r'(^LOC|LINC)[1-9]*')])
remove_genes = ncRNA_genes + LINC_genes
genes_keep = [gene for gene in rna.var_names if gene not in remove_genes]
rna = rna[:,genes_keep].copy()

rna.obs['Age_group'] = rna.obs['Age'].apply(lambda x: 'Young' if x <= 50 else 'Old')
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)


# Figure 3K
G = create_eregulon_network(
    cell_type='CD8_CTL_GZMK',
    cell_type_col='cell_type_l4',
    tf_set={'NFATC3', 'E2F3', 'SMAD3', 'MAF', 'RORA', 
           'NFATC2', 'EOMES', 'GFI1', 'SP4', 'MYBL1'},
    rna_data=rna,
    peak_df=peak_df,
    eRegulon_metadata=eRegulon_metadata,
    output_path='CD8_CTL_GZMK_network.pdf'
)

# Figure 3L
G = create_eregulon_network(
    cell_type='intMono_GFRA2',
    cell_type_col='cell_type_l4',
    tf_set={'KLF11','PLAGL2','CUX1','ZBTB7A','MBD2'},
    rna_data=rna,
    peak_df=peak_df,
    eRegulon_metadata=eRegulon_metadata,
    output_path='intMono_GFRA2_network.pdf'
)
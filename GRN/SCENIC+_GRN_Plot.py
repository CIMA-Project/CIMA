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
    

    cell_type_data = rna_data[rna_data.obs[cell_type_col] == cell_type].copy()
    sc.tl.rank_genes_groups(cell_type_data, group_col, method="wilcoxon")
    DEG_df = sc.get.rank_genes_groups_df(cell_type_data,None)

    DEG_sig_df = DEG_df[(DEG_df['group'] == group_value) & 
                       (DEG_df['pvals_adj'] < deg_pval_thresh) & 
                       (DEG_df['logfoldchanges'] > deg_logfc_thresh)]
    DEG_sig_list = set(DEG_sig_df['names'])
    

    peak_set = set(peak_df[peak_df[cell_type] == True].index)
    

    filtered_eRegulon = eRegulon_metadata[
        (eRegulon_metadata[tf_col].isin(tf_set)) &
        (eRegulon_metadata[gene_col].isin(DEG_sig_list)) &
        (eRegulon_metadata[region_col].isin(peak_set))
    ].copy()
    

    filtered_eRegulon[tf_col] = tf_prefix + filtered_eRegulon[tf_col]
    

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
    

    node_shapes = {'TF': 'h', 'Region': 'D', 'Gene': 'o'}
    node_colors = {'Region': 'blue'}
    node_sizes = {'TF': 1800, 'Region': 60, 'Gene': 300}

    node_type = nx.get_node_attributes(G, 'type')
    node_labels = nx.get_node_attributes(G, 'label')
    

    tf_labels = list(set(nx.get_edge_attributes(G, 'tf').values()))
    color_map = plt.cm.get_cmap('tab20', len(tf_labels))
    tf_color_map = {tf: color_map(i) for i, tf in enumerate(tf_labels)}
    

    def layered_layout(G, layers):
        pos = {}
        for i, (layer, nodes) in enumerate(layers.items()):
            angle_step = 2 * np.pi / len(nodes)
            radius = (i + 1) * 1.5
            for j, node in enumerate(nodes):
                angle = j * angle_step
                pos[node] = (radius * np.cos(angle), radius * np.sin(angle))
        return pos
    

    node_type = nx.get_node_attributes(G, 'type')
    layers = {
        'TF': [n for n in G.nodes if node_type[n] == 'TF'],
        'Region': [n for n in G.nodes if node_type[n] == 'Region'],
        'Gene': [n for n in G.nodes if node_type[n] == 'Gene']
    }
    

    pos = layered_layout(G, layers)
    

    vmin = DEG_sig_df['logfoldchanges'].min()
    vmax = DEG_sig_df['logfoldchanges'].max()
    norm = mcolors.TwoSlopeNorm(vcenter=0, vmax=vmax)
    cmap = plt.cm.bwr
    
    gene_colors = []
    for gene in layers['Gene']:
        gene_name = gene
        if gene_name in deg_df['names'].values:
            logfc = deg_df[deg_df['names'] == gene_name]['logfoldchanges'].values[0]
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
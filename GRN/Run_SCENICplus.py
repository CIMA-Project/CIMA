import dill
import scanpy as sc
import os
import warnings
warnings.filterwarnings("ignore")
import pandas
import pyranges
# Set stderr to null to avoid strange messages from ray
import sys
import ray
ray._private.ray_constants.DEFAULT_OBJECT_STORE_MEMORY_PROPORTION = 0.85
ray._private.ray_constants.DEFAULT_OBJECT_STORE_MAX_MEMORY_BYTES = 1.4*10^12
_stderr = sys.stderr
null = open(os.devnull,'wb')
work_dir = 'CIMA/scenicplus/'
tmp_dir = 'CIMA/Temp'

rna = sc.read_h5ad('CIMA/Pseudo_multiomics/CIMA_scRNA_Pseudo_multiomics.h5ad')
rna.obs.index = rna.obs['new_barcode']
cistopic_obj = dill.load(open(os.path.join(work_dir, 'pycisTopic/CIMA_scATAC_PeakRecalling_Metacell_cisTopic_obj.pkl'), 'rb'))

ncRNA_genes = list(rna.var_names[rna.var_names.str.match(r'^[A-Z][A-Z][0-9].*\.[0-9]')])
LINC_genes = list(rna.var_names[rna.var_names.str.match(r'(^LOC|LINC)[1-9]*')])
remove_genes = ncRNA_genes + LINC_genes
genes_keep = [gene for gene in rna.var_names if gene not in remove_genes]
len(genes_keep)
rna = rna[:,genes_keep].copy()

def list_pkl_files(directory):
    
    files_and_dirs = os.listdir(directory)
    
    pkl_files = [f for f in files_and_dirs if f.endswith('.pkl')]
    
    return pkl_files

# Load cistarget and DEM motif enrichment results
motif_enrichment_dict={}
import pickle
from pycistarget.motif_enrichment_dem import *
path = os.path.join(work_dir, 'motifs_custom')
pkl_files = list_pkl_files(path)

for pkl in pkl_files:
    infile = open(path + pkl, 'rb')
    name = os.path.splitext(pkl)[0]
    motif_enrichment_dict[name] = pickle.load(infile)
    infile.close()

from scenicplus.scenicplus_class import create_SCENICPLUS_object
import numpy as np
scplus_obj = create_SCENICPLUS_object(
    GEX_anndata = rna,
    cisTopic_obj = cistopic_obj,
    menr = motif_enrichment_dict,
    ACC_prefix = 'ACC_',
    GEX_prefix = 'GEX_',
    bc_transform_func = lambda x: x,
    normalize_imputed_acc = False
)

ensembl_version_dict = {'105': 'http://www.ensembl.org',
                        '104': 'http://may2021.archive.ensembl.org/',
                        '103': 'http://feb2021.archive.ensembl.org/',
                        '102': 'http://nov2020.archive.ensembl.org/',
                        '101': 'http://aug2020.archive.ensembl.org/',
                        '100': 'http://apr2020.archive.ensembl.org/',
                        '99': 'http://jan2020.archive.ensembl.org/',
                        '98': 'http://sep2019.archive.ensembl.org/',
                        '97': 'http://jul2019.archive.ensembl.org/',
                        '96': 'http://apr2019.archive.ensembl.org/',
                        '95': 'http://jan2019.archive.ensembl.org/',
                        '94': 'http://oct2018.archive.ensembl.org/',
                        '93': 'http://jul2018.archive.ensembl.org/',
                        '92': 'http://apr2018.archive.ensembl.org/',
                        '91': 'http://dec2017.archive.ensembl.org/',
                        '90': 'http://aug2017.archive.ensembl.org/',
                        '89': 'http://may2017.archive.ensembl.org/',
                        '88': 'http://mar2017.archive.ensembl.org/',
                        '87': 'http://dec2016.archive.ensembl.org/',
                        '86': 'http://oct2016.archive.ensembl.org/',
                        '80': 'http://may2015.archive.ensembl.org/',
                        '77': 'http://oct2014.archive.ensembl.org/',
                        '75': 'http://feb2014.archive.ensembl.org/',
                        '54': 'http://may2009.archive.ensembl.org/'}

import pybiomart as pbm
def test_ensembl_host(scplus_obj, host, species):
    dataset = pbm.Dataset(name=species+'_gene_ensembl',  host=host)
    annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
    annot.columns = ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    annot['Chromosome'] = annot['Chromosome'].astype('str')
    filter = annot['Chromosome'].str.contains('CHR|GL|JH|MT')
    annot = annot[~filter]
    annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    gene_names_release = set(annot['Gene'].tolist())
    ov=len([x for x in scplus_obj.gene_names if x in gene_names_release])
    print('Genes recovered: ' + str(ov) + ' out of ' + str(len(scplus_obj.gene_names)))
    return ov

n_overlap = {}
for version in ensembl_version_dict.keys():
    print(f'host: {version}')
    try:
        n_overlap[version] =  test_ensembl_host(scplus_obj, ensembl_version_dict[version], 'hsapiens')
    except:
        print('Host not reachable')
v = sorted(n_overlap.items(), key=lambda item: item[1], reverse=True)[0][0]
print(f"version: {v} has the largest overlap, use {ensembl_version_dict[v]} as biomart host")

biomart_host = "http://sep2019.archive.ensembl.org/"

scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())

from scenicplus.wrappers.run_scenicplus import run_scenicplus
try:
    sys.stderr = open(os.devnull, "w")  # silence stderr
    run_scenicplus(
        scplus_obj = scplus_obj,
        variable = ['GEX_cell_type_l4'],
        species = 'hsapiens',
        assembly = 'hg38',
        tf_file = 'utoronto_human_tfs_v_1.01.txt',
        save_path = os.path.join(work_dir, 'scenicplus_custom_cell_type_l4'),
        biomart_host = biomart_host,
        upstream = [1000, 150000],
        downstream = [1000, 150000],
        calculate_TF_eGRN_correlation = False,
        calculate_DEGs_DARs = False,
        export_to_loom_file = False,
        export_to_UCSC_file = False,
        path_bedToBigBed = './',
        n_cpu = 96,
        _temp_dir = os.path.join(tmp_dir, 'ray_spill'),
        **{ '_system_config':{
        # Allow spilling until the local disk is 99% utilized.
        # This only affects spilling to the local file system.
        "local_fs_capacity_threshold": 0.99},"dashboard_host":"0.0.0.0"})
    sys.stderr = sys.__stderr__  # unsilence stderr
except Exception as e:
    #in case of failure, still save the object
    dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus_pseudo_multiomics_obj.pkl'), 'wb'), protocol=-1)
    raise(e)

dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus_pseudo_multiomics_obj.pkl'), 'wb'), protocol=-1)
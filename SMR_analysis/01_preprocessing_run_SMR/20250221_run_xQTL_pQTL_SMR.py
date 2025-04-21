import argparse
import os 
import pandas as pd

parser = argparse.ArgumentParser(description='Define arguments')
parser.add_argument('--CT', type=str, required=True)

args = parser.parse_args()
print(args.CT)

celltype = args.CT

eQTL_CT = pd.read_csv('/CIMA/Data/20241230_xQTL_cell_sample_number/20241230_celltype_70_scRNA.csv')
eQTL_CT = eQTL_CT['final_annotation'].to_list()
caQTL_CT = pd.read_csv('/CIMA/Data/20241230_xQTL_cell_sample_number/20241230_celltype_70_scATAC.csv')
caQTL_CT = caQTL_CT['final_annotation'].to_list()
overlap_CT = list(set(eQTL_CT).intersection(caQTL_CT))

#opera = '/media/scPBMC1_AnalysisDisk1/huangzhuoli/Script_HPC/genetics_software/opera_Linux'
smr = '/media/scPBMC1_AnalysisDisk1/huangzhuoli/Script_HPC/software_gaoyue/SMR/smr-1.3.1-linux-x86_64/smr-1.3.1'

#OPERA_output_dir_1 = '/CIMA/Result/xQTL_GWAS_OPERA/stage_1/BBJ/'
#OPERA_output_dir_2 = '/CIMA/Result/xQTL_GWAS_OPERA/stage_2/BBJ/'
SMR_output_dir = '/CIMA/Result/xQTL_GWAS_SMR/pQTL/'
#xQTL_list_dir = '/CIMA/Data/xQTL_list_for_opera/'

GWAS_dir = '/CIMA/Data/pQTL_ma/'
eQTL_output_dir = '/CIMA/Data/SMR_eQTL_besd/besd/'
caQTL_output_dir = '/CIMA/Data/SMR_caQTL_besd/besd/'
eQTL = f'{eQTL_output_dir}{celltype}'
caQTL = f'{caQTL_output_dir}{celltype}'
bim = '/CIMA/genetics/qc/10.maf01'

trait_all = pd.read_csv('/CIMA/Data/pQTL_ma/pQTL_list.txt',sep='\t',header=None)
trait_all = trait_all[0].to_list()

for trait in trait_all:
    print(trait)
    #eQTL
    if not os.path.exists(f'{SMR_output_dir}{trait}'):
        # 如果文件夹不存在，创建文件夹
        os.makedirs(f'{SMR_output_dir}{trait}')

    if celltype in eQTL_CT :
        gi = f'{smr} --bfile {bim} --gwas-summary {GWAS_dir}{trait}.ma --beqtl-summary {eQTL} --thread-num 1 --diff-freq-prop 0.5 --out {SMR_output_dir}{trait}/eQTLsmr_{celltype}'
        os.system(gi)

    #caQTL
    if celltype in caQTL_CT :
        gi = f'{smr} --bfile {bim} --gwas-summary {GWAS_dir}{trait}.ma --beqtl-summary {caQTL}  --thread-num 1 --diff-freq-prop 0.5 --out {SMR_output_dir}{trait}/caQTLsmr_{celltype}'
        os.system(gi)

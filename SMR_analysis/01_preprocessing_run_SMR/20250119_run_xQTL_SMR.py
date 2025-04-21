import argparse
import os 

parser = argparse.ArgumentParser(description='Define arguments')
parser.add_argument('--CT', type=str, required=True)

args = parser.parse_args()
print(args.CT)

celltype = args.CT

eQTL_output_dir = '/CIMA/Data/SMR_eQTL_besd/besd/'
caQTL_output_dir = '/CIMA/Data/SMR_caQTL_besd/besd/'
SMR_output_dir = '/CIMA/Result/xQTL_SMR/'

smr = '/media/scPBMC1_AnalysisDisk1/huangzhuoli/Script_HPC/software_gaoyue/SMR/smr-1.3.1-linux-x86_64/smr-1.3.1'
bim = '/CIMA/genetics/qc/10.maf01'

eQTL = f'{eQTL_output_dir}{celltype}'
caQTL = f'{caQTL_output_dir}{celltype}'

gi = f'{smr} --bfile {bim} --beqtl-summary {caQTL} --beqtl-summary {eQTL} --thread-num 1 --out {SMR_output_dir}{celltype}'
os.system(gi)


import argparse
import os 
import pandas as pd

parser = argparse.ArgumentParser(description='Define arguments')
parser.add_argument('--CT', type=str, required=True)

args = parser.parse_args()
print(args.CT)

celltype = args.CT

smr = '/media/scPBMC1_AnalysisDisk1/huangzhuoli/Script_HPC/software_gaoyue/SMR/smr-1.3.1-linux-x86_64/smr-1.3.1'
SMR_output_dir = '/CIMA/Result/dynamic/SMR_BBJ/'


GWAS_dir = '/CIMA/Data/BBJ_GWAS_for_SMR_OPERA/ma_file/'
eQTL_output_dir = '/CIMA/Data/dynamic/for_SMR/besd/'
eQTL = f'{eQTL_output_dir}{celltype}'
bim = '/CIMA/genetics/qc/10.maf01'

trait_all = pd.read_csv('/CIMA/Data/BBJ_GWAS_for_SMR_OPERA/BBJ_disease.txt',sep='\t',header=None)
trait_all = trait_all[0].to_list()

for trait in trait_all:
    print(trait)
    #eQTL
    if not os.path.exists(f'{SMR_output_dir}{trait}'):
        # 如果文件夹不存在，创建文件夹
        os.makedirs(f'{SMR_output_dir}{trait}')
    gi = f'{smr} --bfile {bim} --gwas-summary {GWAS_dir}hum0197.v3.BBJ.{trait}.v1.ma --beqtl-summary {eQTL} --thread-num 1  --out {SMR_output_dir}{trait}/eQTLsmr_{celltype}'
    os.system(gi)
#!/bin/bash

# 使用cut和tail命令提取第一列并跳过第一行
result=$(cut -d ',' -f 1 //CIMA/Data/20241230_xQTL_cell_sample_number/20241230_celltype_70_scRNA.csv | tail -n +2)

# 逐行处理结果
while IFS= read -r line; do
    # 打印当前时间和每一行的值
    echo "$(date) - Running script for cell: $line"
    
    # 运行python脚本
    /home/huangzhuoli/mambaforge/envs/tensorqtl/bin/python 20250103_run_cis_eQTL.py --cell "$line"
    
    echo "$(date) - Finished script for cell: $line"
done <<< "$result"

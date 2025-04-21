#!/bin/bash

# 定义需要的参数列表
cell_types=('CD4_T' 'CD8_T' 'NK' 'Myeloid' 'B')

# 逐个处理这些值
for cell_type in "${cell_types[@]}"; do
    # 打印当前时间和每一行的值
    echo "$(date) - Running script for cell: $cell_type"

    # 运行python脚本
    /home/huangzhuoli/mambaforge/envs/tensorqtl/bin/python 20250214_run_cis_eQTL_L1.py --cell "$cell_type"

    echo "$(date) - Finished script for cell: $cell_type"
done

#!/bin/bash

# 定义需要的参数列表
cell_types=('CD4_T' 'CD8_T' 'NK' 'Myeloid' 'B')

# 逐个处理这些值
for cell_type in "${cell_types[@]}"; do
    # 打印每个值
    echo "$cell_type"

    # 执行Rscript命令并将输出重定向到相应的文件，使用nohup使命令在后台运行
    nohup /home/huangzhuoli/mambaforge/envs/peer/bin/Rscript ./20250214_cal_peer_eQTL.R "$cell_type" > "/CIMA/Data/eQTL_L1/peer/peer_out/$cell_type.out" &
done

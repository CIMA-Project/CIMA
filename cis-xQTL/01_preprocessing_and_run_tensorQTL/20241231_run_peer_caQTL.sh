#!/bin/bash

# 使用cut和tail命令提取第一列并跳过第一行
result=$(cut -d ',' -f 1 /CIMA/Data/20241230_xQTL_cell_sample_number/20241230_celltype_70_scATAC.csv | tail -n +2)

# 逐行处理结果
while IFS= read -r line; do
    # 打印每一行的值
    echo "$line"
    
    # 执行Rscript命令并将输出重定向到相应的文件，使用nohup使命令在后台运行
    nohup /home/huangzhuoli/mambaforge/envs/peer/bin/Rscript ./20241231_cal_peer_caQTL.R "$line" > "/CIMA/Data/caQTL/peer/peer_out/$line.out" &
done <<< "$result"

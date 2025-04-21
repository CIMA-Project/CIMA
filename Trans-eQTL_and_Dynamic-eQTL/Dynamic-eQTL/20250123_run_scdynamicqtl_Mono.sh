#!/bin/bash
cat Mono_cell_dynamic_eGene.txt | xargs -I {} -P 60 /home/zhengyuhui/mambaforge/envs/pbmc/lib/R/bin/Rscript 20250123_run_scdynamicqtl.r {}

#!/bin/bash
seq 1 321 | xargs -n 1 -P 70 bash ./20250114_meta_gwas.sh
seq 1 653 | xargs -n 1 -P 70 bash ./20250114_lipid_gwas.sh

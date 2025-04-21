#!/bin/bash
cat rest_6CT.txt | xargs -I {} -P 60 python3 20250221_run_xQTL_pQTL_SMR.py --CT {}

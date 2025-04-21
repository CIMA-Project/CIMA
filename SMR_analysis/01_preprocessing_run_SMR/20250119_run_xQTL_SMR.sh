#!/bin/bash
cat overlap_42_CT.txt | xargs -I {} -P 50 python3 20250119_run_xQTL_SMR.py --CT {}

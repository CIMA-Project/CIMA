#!/bin/bash
cat all_69_CT.txt | xargs -I {} -P 70 python3 20250119_run_xQTL_lime_SMR.py --CT {}

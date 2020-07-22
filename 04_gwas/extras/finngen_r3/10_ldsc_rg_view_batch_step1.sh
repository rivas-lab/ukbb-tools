#!/bin/bash
set -beEuo pipefail

log_dir="/scratch/groups/mrivas/public_data/summary_stats/finngen_r3/ldsc/UKB_WB_rg"
out_dir="/scratch/groups/mrivas/public_data/summary_stats/finngen_r3/ldsc/UKB_WB_rg.20200722-093650"
ml load ukbb-tools

find ${log_dir} -name "*.log" | sort --parallel 6 -k1,1V -k2,2V > ${out_dir}/LDSC.rg.lst

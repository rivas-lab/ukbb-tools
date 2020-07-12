#!/bin/bash
set -beEuo pipefail

log_dir="/scratch/groups/mrivas/public_data/summary_stats/finngen_r3/ldsc/UKB_WB_rg"
out_f="/oak/stanford/groups/mrivas/public_data/summary_stats/finngen_r3/UKB_WB_rg.$(date +%Y%m%d-%H%M%S).tsv"

ml load ukbb-tools

find ${log_dir} -name "*.log" | sort --parallel 6 -k1,1V -k2,2V | ldsc_rg_view.sh -l /dev/stdin \
| sed -e "s%/scratch/groups/mrivas/public_data/summary_stats/finngen_r3/summary_stats_hg19_ldsc_munge/finngen_r3_%%g" \
| sed -e "s/.hg19.sumstats.gz//g" \
| sed -e "s%/oak/stanford/groups/mrivas/ukbb24983/array-combined/ldsc/white_british/ukb24983_v2_hg19.%%g" \
| sed -e "s/.array-combined.sumstats.gz//g" > ${out_f}

echo ${out_f}

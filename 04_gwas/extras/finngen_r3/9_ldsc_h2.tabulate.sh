#!/bin/bash
set -beEuo pipefail

in_d="/scratch/groups/mrivas/public_data/summary_stats/finngen_r3/ldsc/h2"
out_f="9_ldsc_h2.tsv"

find ${in_d} -type f -name "*.log" | sort | bash ../../../07_LDSC/helpers/ldsc_h2_view.sh -l /dev/stdin | sed -e "s/.hg19.sumstats.gz//g" | sed -e "s%/scratch/groups/mrivas/public_data/summary_stats/finngen_r3/summary_stats_hg19_ldsc_munge/finngen_r3_%%g" | tee ${out_f}

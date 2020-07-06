#!/bin/bash
set -beEuo pipefail

in_d="/oak/stanford/groups/mrivas/ukbb24983/array-combined/ldsc/h2"
out_f="2_ldsc_h2.white_british.tsv"

find ${in_d} -type f -name "*.log" | sort | bash ../../../07_LDSC/helpers/ldsc_h2_view.sh -l /dev/stdin | sed -e "s/.array-combined.sumstats.gz//g" | sed -e "s%/oak/stanford/groups/mrivas/ukbb24983/array-combined/ldsc/white_british/ukb24983_v2_hg19.%%g" | tee ${out_f}


#!/bin/bash
set -beEuo pipefail

pop=$1
in_d="/oak/stanford/groups/mrivas/ukbb24983/array-combined/ldsc/h2"
out_f="2_ldsc_h2.${pop}.tsv"

find ${in_d} -type f -name "${pop}*.log" | sort | bash ../../../07_LDSC/helpers/ldsc_h2_view.sh -l /dev/stdin | sed -e "s/.array-combined.sumstats.gz//g" | sed -e "s%/oak/stanford/groups/mrivas/ukbb24983/array-combined/ldsc/${pop}/ukb24983_v2_hg19.%%g" | tee ${out_f}

echo ${out_f}

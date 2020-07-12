#!/bin/bash
set -beEuo pipefail

out_dir="/scratch/groups/mrivas/public_data/summary_stats/finngen_r3/ldsc/UKB_WB_rg.20200711-171504"
n_batchs=1000
if [ $# -gt 0 ] ; then n_batchs=$1 ; fi

out_f="${out_dir}/LDSC.rg.tsv"

{
    ! cat "${out_dir}/LDSC.rg.batch1.tsv" | egrep '^#'
    ! seq ${n_batchs} | while read batch_idx ; do cat "${out_dir}/LDSC.rg.batch${batch_idx}.tsv" | egrep -v '^#' ; done
} > ${out_f}

echo ${out_f}

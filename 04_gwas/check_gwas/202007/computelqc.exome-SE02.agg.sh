#!/bin/bash
set -beEuo pipefail

ml load parallel
data_d='/scratch/groups/mrivas/ukbb24983/exome/gwas-qc-SE02'

pop=$1

out_tbl="${data_d}/${pop}.qc.tsv"

{
    echo "#GBE_ID population freq_bin lambda_gc" | tr ' ' '\t'
    
!   find ${data_d}/${pop} -name "${pop}.*.qc.txt" | sort -V | parallel --eta -k "cat {} | grep -v '^#'"
} > ${out_tbl}

echo ${out_tbl}

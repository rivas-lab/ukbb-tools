#!/bin/bash
set -beEuo pipefail

data_d='/scratch/groups/mrivas/ukbb24983/array-combined/gwas-qc-SE02'

pop=$1

out_tbl="${data_d}/${pop}.cnt.tsv"

{
    echo "#GBE_ID population n_lines n_non_NA_lines n_hits n_ld_indep_hits" | tr ' ' '\t'
    
!   find ${data_d}/${pop} -name "${pop}.*.cnt.tsv" | sort -V | parallel --eta -k "cat {} | grep -v '^#'"
} > ${out_tbl}

echo ${out_tbl}

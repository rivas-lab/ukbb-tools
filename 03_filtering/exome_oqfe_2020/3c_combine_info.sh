#!/bin/bash
set -beEuo pipefail

data_d="/scratch/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/download"
idx_f=${data_d}/pvcf_blocks.txt

{
    zstdcat $(dirname ${data_d})/pvcf_pvar/ukb23156_cY_b0_v1.pvar.zst | egrep -v '^##' \
        | egrep '#' | sed -e 's/ID/pvcf_ID/g' | cut -f1-7
    for idx in $(seq 1 977) ; do

        c=$(cat ${idx_f} | awk -v idx=${idx} '($1 == idx){print $2}' | sed -e 's/23/X/g' | sed -e 's/24/Y/g')
        b=$(cat ${idx_f} | awk -v idx=${idx} '($1 == idx){print $3}')
        pvar="$(dirname ${data_d})/pvcf_pvar/ukb23156_c${c}_b${b}_v1.pvar.zst"

        zstdcat ${pvar} | egrep -v '^#' | cut -f1-7
    done

} > $(dirname ${data_d})/ukb23156_pvcf_info.tsv

bgzip -@6 -l9 $(dirname ${data_d})/ukb23156_pvcf_info.tsv

# later, this file is copied to /oak/stanford/groups/mrivas/ukbb24983/exome/qc/oqfe_2020/intermediate_files/ukb23156_pvcf_info.tsv.gz

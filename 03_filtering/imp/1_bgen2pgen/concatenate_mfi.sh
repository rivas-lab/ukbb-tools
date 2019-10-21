#!/bin/bash
set -beEuo pipefail

mfi_dir="/oak/stanford/groups/mrivas/ukbb24983/imp/mfi"

{
    echo "#ID UKB_VAR_ID ORIGINAL_VAR_ID AF_A1 INFO" | tr " " "\t"

    for c in $(seq 1 22) X XY ; do
        zstdcat ${mfi_dir}/ukb_mfi_chr${c}_v3.txt.zst \
        | awk -v FS='\t' -v OFS='\t' -v sep=':' -v c="${c}" \
        '{if($7==$5){af=$6}else{af=1-$6} print c sep $3 sep $4 sep $5, $1, $2, af, $8}'
    done
} | zstd -15 > ${mfi_dir}/ukb_mfi_v3.tsv.zst


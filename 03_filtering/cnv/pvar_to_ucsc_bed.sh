#!/bin/bash
set -beEuo pipefail

cat_or_zcat () {
    local file=$1
    if [ "${file%.gz}.gz" == "${file}" ] || [ "${file%.bgz}.bgz" == "${file}" ] ; then
        zcat ${file}
    elif [ "${file%.zst}.zst" == "${file}" ] ; then
        zstdcat ${file}
    else
        cat ${file}
    fi
}

cnv_pvar="/oak/stanford/groups/mrivas/ukbb24983/cnv/pgen/cnv.pvar"
cnv_ucsc_bed="$(dirname ${cnv_pvar})/cnv.ucsc.bed"

cat_or_zcat ${cnv_pvar} | awk -v FS='\t' '(NR>1){print $3}' | tr ':' '\t' | tr '-' '\t' | tr '_' '\t' | awk -v FS='\t' -v OFS='\t' '{print $1, $2, $3}' | paste /dev/stdin <(cat_or_zcat ${cnv_pvar} | awk -v FS='\t' '(NR>1){print $3}') > ${cnv_ucsc_bed}

echo ${cnv_ucsc_bed}

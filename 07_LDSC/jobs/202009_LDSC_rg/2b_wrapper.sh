#!/bin/bash
set -beEuo pipefail

idx=$1
if [ $# -gt 1 ] ; then offset=$2 ; else offset=0 ; fi
pair_idx_f="GBE_ID.pairs.tsv.gz"

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

line=$(perl -e "print($idx + $offset)")

GBE_ID1=$(cat_or_zcat ${pair_idx_f} | awk -v l=$line '(NR == l){print $1}')
GBE_ID2=$(cat_or_zcat ${pair_idx_f} | awk -v l=$line '(NR == l){print $2}')

bash 1_ldsc_rg.sh ${GBE_ID1} ${GBE_ID2}


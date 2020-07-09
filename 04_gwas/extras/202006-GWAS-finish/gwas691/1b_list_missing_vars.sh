#!/bin/bash
set -beEuo pipefail

pvar="/oak/stanford/groups/mrivas/ukbb24983/array-combined/pgen/ukb24983_cal_hla_cnv.pvar.zst"
vars=$1
out="missing.lst"

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

cat_or_zcat ${pvar} | cut -f3 | sort \
| comm -3 /dev/stdin <(cat_or_zcat ${vars} | cut -f3 | sort) > ${out}

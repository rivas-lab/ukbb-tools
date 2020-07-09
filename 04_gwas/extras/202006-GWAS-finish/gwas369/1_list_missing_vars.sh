#!/bin/bash
set -beEuo pipefail

pvar="/oak/stanford/groups/mrivas/ukbb24983/array-combined/pgen/ukb24983_cal_hla_cnv.pvar.zst"
sumstats="/oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/2007183/40831/others/ukb24983_v2_hg19.INI22149.array-combined.glm.linear.gz"
out="missing.369.lst"

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
| comm -3 /dev/stdin <(cat_or_zcat ${sumstats} | cut -f3 | sort) > ${out}

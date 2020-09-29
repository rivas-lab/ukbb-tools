#!/bin/bash
set -beEuo pipefail

ldmap_d=/scratch/groups/mrivas/ukbb24983/array_imp_combined_no_cnv/ldmap/ldmap_20200928

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

ml load zstd

cat merge.lst.v2.imp.wo.HWE.filter.tsv | cut -f1 | while read pvar ; do
    cat_or_zcat ${pvar} | egrep -v '#' | awk '{print $3}'
done > ${ldmap_d}/ukb24983_cal_hla_imp.variants.lst

echo ${ldmap_d}/ukb24983_cal_hla_imp.variants.lst

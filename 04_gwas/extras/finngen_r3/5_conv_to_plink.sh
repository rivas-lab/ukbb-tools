#!/bin/bash
set -beEuo pipefail

in_f=$1
out_f=$(dirname $(dirname $in_f))/summary_stats_hg19_plink/$(basename ${in_f%.hg19.gz}.hg19.plink.gz)

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

if [ ! -f ${out_f}.tbi ] ; then

    {
        # header
        echo "#CHROM POS ID REF ALT A1 FIRTH? TEST OBS_CT OR LOG(OR)_SE Z_STAT P ERRCODE" | tr " " "\t"

        # body
        cat_or_zcat $in_f | egrep -v '^#' | awk -v FS='\t' -v OFS='\t' '{ print $1, $2, $5, $3, $4, $4, "N", "ADD", 135638, $7, $8, log($7)/$8, $6, "."}'

    } | bgzip -l9 -@2 > ${out_f}

    tabix -f -c '#' -s 1 -b2 -e2 ${out_f}

    echo ${out_f}
fi

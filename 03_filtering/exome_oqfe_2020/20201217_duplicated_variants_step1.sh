#!/bin/bash
set -beEuo pipefail

in_pvar="/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/ukb24983_exomeOQFE.pvar.zst"
out_tsv="$(dirname ${in_pvar})/ukb24983_exomeOQFE.duplicates.noID.tsv"
# later, we moved the files to /oak/stanford/groups/mrivas/ukbb24983/exome/qc/oqfe_2020/intermediate_files/ukb24983_exomeOQFE.duplicates.tsv.gz

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

echo "#CHROM POS REF ALT n" | tr ' ' '\t' > ${out_tsv}

cat_or_zcat ${in_pvar} \
| awk -v OFS='_' '(NR>1){print $1, $2, $4, $5}' \
| sort --parallel 6 | uniq -c \
| awk -v OFS='\t' '($1>1){print $2, $1}' | tr '_' '\t' >> ${out_tsv}

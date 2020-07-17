#!/bin/bash
set -beEuo pipefail

check_file="2_gwas.jobs.status.20200716-232605.tsv"

cat ${check_file} |
awk -v FS='\t' -v OFS='\t' -v T="TRUE" 'NR>1 && ($5 == T && $6 == T){print $4, $7, $8}' |
while read src dst dst_syml ; do

    if [ ! -d $(dirname ${dst}) ] ; then
        mkdir -p $(dirname ${dst})
    fi

    mv ${src} ${dst}
    ln -s ${dst} ${dst_syml}

done

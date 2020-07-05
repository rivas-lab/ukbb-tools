#!/bin/bash
set -beEuo pipefail

in_f="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/extras/202006-GWAS-finish/gwas-current-gz-wc.20200704-155715.annotated.tsv"

if [ $# -gt 0 ] ; then
    in_f=$1
fi

cat $in_f | awk -v FS='\t' -v OFS='\t' '(NR > 1 || $NF == 1080969){ print $1, $2, $3 }' | while read GBE_ID pop sumstats_symlink ; do

    LDSC_f=/oak/stanford/groups/mrivas/ukbb24983/array-combined/ldsc/${pop}/ukb24983_v2_hg19.${GBE_ID}.array-combined.sumstats.gz

    if [ ! -f ${LDSC_f} ] ; then
        echo $sumstats_symlink
    fi
done

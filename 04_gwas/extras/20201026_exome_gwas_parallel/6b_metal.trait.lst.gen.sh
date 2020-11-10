#!/bin/bash
set -beEuo pipefail

out_f="6b_metal.trait.lst.gen.$(date +%Y%m%d-%H%M%S).lst"

cat /scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l.20201105-081350.tsv | awk '(NR>1){print $2}' | sort -uV \
    | while read GBE_ID ; do
    f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/metal/ukb24983_exomeOQFE.${GBE_ID}
    flog=${f}.metal.info.txt
    fgwas=${f}.metal.tsv.gz
    
    if [ ! -s ${flog} ] && [ ! -s ${fgwas} ] ; then
        echo $GBE_ID
    fi
done > ${out_f}

echo ${out_f}


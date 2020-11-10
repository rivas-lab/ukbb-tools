#!/bin/bash
set -beEuo pipefail

out_f="6b_metal.trait.lst.gen.$(date +%Y%m%d-%H%M%S).lst"

bash 5b_finished_list.sh /dev/stdout | egrep -v '#' | cut -f2 | sort | uniq -c | awk '($1==7){print $2}' \
    | while read GBE_ID ; do
    f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/metal/ukb24983_exomeOQFE.${GBE_ID}
    flog=${f}.metal.info.txt
    fgwas=${f}.metal.tsv.gz
    
    if [ ! -s ${flog} ] && [ ! -s ${fgwas} ] ; then
        echo $GBE_ID
    fi
done > ${out_f}

echo ${out_f}


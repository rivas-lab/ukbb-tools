#!/bin/bash
set -beEuo pipefail

#wc_l_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l.20201104-133054.tsv
wc_l_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l.20201104-162723.tsv

copy_for_bkup () {
    local f=$1
    
    bk_f=$(dirname $(dirname $f))/patch_fix_bkup/$(basename $(dirname $f))/$(basename $f) 

    if [ ! -d $(dirname $bk_f) ] ; then mkdir -p $(dirname $bk_f) ; fi
    
#    if [ ! -s $bk_f ] ; then
        cp $f $bk_f
#    fi
}

export -f copy_for_bkup

cat ${wc_l_f} | awk -v FS='\t' '( (NR>1) && (! ($4 == 1 && ($5 == 17777817 || $5 == 17777951)) ) ){print $NF}' \
| parallel -j6 --eta 'copy_for_bkup {}'

# cat ${wc_l_f} | awk -v FS='\t' '(NR>1&&$4 != 1){print $NF}' \
# | parallel -j6 --eta 'copy_for_bkup {}'

# cat ${wc_l_f} | awk -v FS='\t' '(NR>1&&$4 != 1){print $NF}' \
# | while read f ; do echo $f ; copy_for_bkup $f ; done

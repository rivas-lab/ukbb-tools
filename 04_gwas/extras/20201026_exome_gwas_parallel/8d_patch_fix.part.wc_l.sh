#!/bin/bash
set -beEuo pipefail

# wc_l_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l.20201104-090826.tsv
# wc_l_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l.20201104-121324.tsv
wc_l_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l.20201104-133054.tsv

# cat ${wc_l_f} | awk -v FS='\t' '(NR>1&&$4 != 1){print $NF}'
cat ${wc_l_f} | awk -v FS='\t' '( (NR>1) && (! ($4 == 1 && ($5 == 17777817 || $5 == 17777951)) ) ){print $1, $2}' \
| while read pop GBE_ID ; do
        find /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/${pop}-batch -name "*.${GBE_ID}.glm.logistic.hybrid*" -o -name "*.${GBE_ID}.glm.linear*"  | sort -V 
done
exit 0
############
bash 8d_patch_fix.part.wc_l.sh > 7b_wc_l_joblist.$(date +%Y%m%d-%H%M%S).lst

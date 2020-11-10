#!/bin/bash
set -beEuo pipefail

# wc_l_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l.20201104-090826.tsv
# wc_l_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l.20201104-121324.tsv
# wc_l_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l.20201104-133054.tsv
# wc_l_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l.20201104-151418.tsv
# wc_l_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l.20201104-162723.tsv
# wc_l_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l.20201104-171349.tsv
wc_l_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l.20201104-175636.tsv

# cat ${wc_l_f} | awk -v FS='\t' '(NR>1&&$4 != 1){print $NF}'
cat ${wc_l_f} | awk -v FS='\t' '( (NR>1) && (! ($4 == 1 && ($5 == 17777817 || $5 == 17777951)) ) ){print $NF}'

exit 0
############
bash 8c_patch_fixed.wc_l.sh > 7b_wc_l_joblist.$(date +%Y%m%d-%H%M%S).lst

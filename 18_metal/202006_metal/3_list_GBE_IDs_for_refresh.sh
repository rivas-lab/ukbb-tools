#!/bin/bash
set -beEuo pipefail

## this script is actually not used

############################################################
# constants
############################################################

f_old="../../04_gwas/extras/202006-GWAS-finish/gwas-current-gz-wc.20200704-155715.annotated.tsv"
f_new="../../04_gwas/extras/202006-GWAS-finish/gwas-current-gz-wc.20200717-173655.annotated.tsv"

out_prefix="3_GBE_IDs.$(date +%Y%m%d-%H%M%S)"

############################################################
# body
############################################################

# find the GBE_IDs to copy
comm -12 <( cat ${f_old} | egrep -v '^#' | tr '\t' ':' | sort ) <( cat ${f_new} | egrep -v '^#' | tr '\t' ':' | sort ) | awk -v FS=':' '{print $1}' | sed -e 's/\s//g' | sort -u > ${out_prefix}.copy_from_old.lst

# find the list of GBE_IDs for re-run.
cat ${f_new} | egrep -v '^#' | cut -f1 | sort -u | comm -23 /dev/stdin ${out_prefix}.copy_from_old.lst > ${out_prefix}.rerun.lst

#
echo "This should produce the empty list"
cat ${f_new} | egrep -v '^#' | cut -f1 | sort -u | comm -3 /dev/stdin <( cat ${out_prefix}.copy_from_old.lst ${out_prefix}.rerun.lst | sort -u)

comm -12 ${out_prefix}.copy_from_old.lst ${out_prefix}.rerun.lst

echo ${out_prefix}.rerun.lst
echo ${out_prefix}.copy_from_old.lst

#!/bin/bash
set -beEuo pipefail

# In the original version of 1_generate_input_list.sh, we incorrectly specified `NR>1 || $NF == 1080969`, but it should have been `NR>1 && $NF == 1080969`. This results resulted in 909 extra munged files.
# Those were NOT used in the heritability analysis. In this script, we remove those 909 files.

cat /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/extras/202006-GWAS-finish/gwas-current-gz-wc.20200704-155715.annotated.tsv | awk -v FS='\t' -v OFS='\t' '(NR > 1 && $NF != 1080969){ print $1, $2, $3 }' | grep -v '_X.array-combined' | while read GBE_ID pop sumstats_symlink ; do

    LDSC_f=/oak/stanford/groups/mrivas/ukbb24983/array-combined/ldsc/${pop}/ukb24983_v2_hg19.${GBE_ID}.array-combined.sumstats.gz

    if [ -f ${LDSC_f} ] ; then
        echo $LDSC_f
    fi
done | tar --remove-files -czvf 1_remove-incomplete-$(date +%Y%m%d-%H%M%S).bkup.tar.gz -T -

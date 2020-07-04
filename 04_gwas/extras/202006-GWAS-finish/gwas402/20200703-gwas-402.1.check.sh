#!/bin/bash
set -beEuo pipefail

# dir402="/oak/stanford/groups/mrivas/users/guhan/sandbox/variants_402"
dir402="/scratch/groups/mrivas/users/ytanigaw/20200703_gwas_merge_bkup/variants_402"

idx_file="gwas-current-gz-wc.20200629-131527.tsv"

# step 1, copy the files to /scratch

# step 2, check the files in the 402 variants dir
# find ${dir402} -type f -name "*.glm.*" | while read f ; do wc $f ; done | awk '$1 != 403'
# this file was empty:
#   white_british/ukb24983_v2_hg19.2401-3562.array-combined.HC976.glm.logistic.hybrid
# I deleted this file

# step 3, check if we have one file per (GBE_ID, pop)
# --> it is not the case
# find ${dir402}/white_british -type f -name "*.glm.*" | while read f ; do basename $f ; done | awk -v FS='.' '{print $2}' | sort | uniq -c
#    1338 1
#    1238 2
#     483 3
#     398 4

# step 4, check if we have results for 402 variants
{
    echo "#GBE_ID population f402"
    find ${dir402} -type f -name "*.glm.*" | while read f ; do
        pop=$(basename $(dirname $f))
        GBE_ID=$(basename $f | awk -v FS='.' '{print $4}')
        echo $GBE_ID $pop $f
    done
} | tr " " "\t" > 20200703-gwas-402.1.check.tsv
exit 0

cat ${idx_file} | awk -v wc_l=1080567 '($3 == wc_l){print $2}' | while read f ; do
    pop=$(basename $(dirname $f))
    GBE_ID=$(basename $f | awk -v FS='.' '{print $2}')
    f402_num=$(find ${dir402}/${pop} -type f -name "*.${GBE_ID}.*.glm.*" | wc -l)
    if [ ${f402_num} -lt 1 ] ; then
        echo $pop $GBE_ID $f
    fi
done | wc

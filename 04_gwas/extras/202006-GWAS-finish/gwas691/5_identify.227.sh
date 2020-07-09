#!/bin/bash
set -beEuo pipefail

{
echo "#GBE_ID pop log"

find data -name "*.log" | sed -e "s%data/%%g" | sed -e "s/.array-combined.log//g" | sort \
| comm -3 /dev/stdin <( find data -name "*glm*" | sed -e "s%data/%%g" | sed -e "s/.array-combined.PHENO1.glm//g" | sed -e "s/.linear$//g" | sed -e "s/.logistic.hybrid$//g" | sort ) \
| sed -e "s%ukb24983_v2_hg19.%%g" | tr '/' '\t' \
| while read pop GBE_ID ; do

echo $GBE_ID $pop data/${pop}/ukb24983_v2_hg19.${GBE_ID}.array-combined.log

done

} | tr ' ' '\t' > 5_identify.227.tsv

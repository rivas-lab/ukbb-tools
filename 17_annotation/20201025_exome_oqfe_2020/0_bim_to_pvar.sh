#!/bin/bash
set -beEuo pipefail

pvar_f=/scratch/groups/mrivas/ukbb24983/exome-oqfe2020-annotation/UKBexomeOQFE.pvar.zst

{
echo "#CHROM POS ID REF ALT" | tr ' ' '\t'
for c in $(seq 1 22) X Y ; do 
cat /scratch/groups/mrivas/ukbb24983/exome-oqfe2020/UKBexomeOQFE_chr${c}.bim \
    | awk -v OFS='\t' '{print $1, $4, $2, $6, $5}'
done
} > ${pvar_f%.zst}

zstd -9 ${pvar_f%.zst}
echo ${pvar_f}


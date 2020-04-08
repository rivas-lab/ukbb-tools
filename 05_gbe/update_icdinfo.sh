#!/bin/bash
set -beEuo pipefail

# please run this script from right here in 05_gbe

. ../04_gwas/04_gwas_misc.sh

# path to extra icdinfo subsets to import (from external sources like MVP, BROAD, etc.)
if [ $# -gt 0 ]; then
  ICD_DIR=$1
else
  ICD_DIR=../04_gwas
fi

# output redirect
if [ $# -gt 1 ]; then
  OUT_DIR=$2
else 
  OUT_DIR=`pwd`
fi

pop="white_british"
if [ $# -gt 2 ] ; then pop=$3 ; fi
field=$(get_field_from_pop $pop)

# updates icdinfo
(
# first take phenotype info
awk -F'\t' -v field=${field} 'BEGIN{OFS="\t"}(NR > 1 && length($1) && length($2) && length($field)){print $1,$field,$2,$field,$field,"Y"}' phenotype_info.tsv ;

# then find additional info
#find ${ICD_DIR} -name "info.tsv" | xargs -i awk 'BEGIN{OFS="\t"}(length($1) && length($2) && length($3)){print $1,$2,$3,$2,$2,"Y"}' {} ; 
) | tr " " "_" | sort -k1,1 -k2,2nr | awk '!_[$1]++' > ${OUT_DIR}/icdinfo.${pop}.txt

if [ "${pop}" == "white_british" ] ; then
    cd ${OUT_DIR}
    ln -sf icdinfo.${pop}.txt icdinfo.txt
fi


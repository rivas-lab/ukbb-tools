#!/bin/bash
set -beEuo pipefail

# please run this script from right here in 05_gbe
. ../04_gwas/04_gwas_misc.sh

# path to extra icdinfo subsets to import (from external sources like MVP, BROAD, etc.)
OUT_DIR=`pwd`

pop="white_british"
# population specification
if [ $# -gt 0 ]; then pop=$1 ; fi
field=$(get_field_from_pop $pop)

if [ -z "$field" ]; then echo "Invalid population!"; exit 1; fi

# updates icdinfo by manipulating phenotype_info.tsv
awk -F'\t' -v field=${field} 'BEGIN{OFS="\t"}(NR > 1 && length($1) && length($2) && length($field)){print $1,$field,$2,$field,$field,"Y"}' phenotype_info.tsv | tr " " "_" | sort -k1,1 -k2,2nr | awk '!_[$1]++' > ${OUT_DIR}/icdinfo.${pop}.txt

if [ "${pop}" == "white_british" ] ; then
    cd ${OUT_DIR}
    ln -sf icdinfo.${pop}.txt icdinfo.txt
fi

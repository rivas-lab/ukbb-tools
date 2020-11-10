#!/bin/bash
set -beEuo pipefail

# constants
data_d=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_filtered_p1e-3

# functions
source 0_functions.sh
export -f cat_or_zcat
export -f master_gwas_dump_file

# cmd args

pop=$1
header_phe=$2

if [ $# -gt 2 ] ; then
    suffix=$3
elif [ "${pop}" == "metal" ] ; then
    suffix="metal.tsv.gz"
else
    suffix=$(get_plink_suffix ${header_phe}).gz
fi

# pop=metal
# header_phe="INI50"
# suffix="metal.tsv.gz"

# show header
find ${data_d}/${pop} -type f -name "*.${header_phe}.*${suffix}" | awk 'NR==1' \
    | while read f ; do cat_or_zcat $f ; done \
    | egrep '^#' | uniq \
    | sed -e 's/^#//g' | sed -e 's/ID/Variant_ID/g' \
    | awk -v OFS='\t' '{print "#population", "GBE_ID", $0}'

# show body
# -- please check the function file (0_functions.sh) for the definition of master_gwas_dump_file
find ${data_d}/${pop} -type f -name "*.${suffix}" | sort -V \
    | parallel -k -j+0 master_gwas_dump_file


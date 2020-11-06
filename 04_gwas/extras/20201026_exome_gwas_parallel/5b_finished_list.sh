#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.0.1"
NUM_POS_ARGS="1"

out_f=${PROGNAME%.sh}.$(date +%Y%m%d-%H%M%S).tsv

if [ $# -gt 0 ] ; then out_f=$1 ; fi

{
echo "#population GBE_ID" | tr ' ' '\t'
for pop in 'white_british' 'non_british_white' 'african' 's_asian' 'e_asian' 'related' 'others' ; do

find /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/${pop} -type f -name "*linear.gz" -o -name "*logistic.hybrid.gz" | awk -v FS='/' '{print $NF}' | awk -v FS='.' '{print $2}' | sort -V \
    | awk -v p=${pop} -v OFS='\t' '{print p, $0}'
done
} > ${out_f}
echo ${out_f}


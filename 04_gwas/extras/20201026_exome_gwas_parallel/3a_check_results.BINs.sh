#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.0.1"
NUM_POS_ARGS="1"

out_f=${PROGNAME%.sh}.$(date +%Y%m%d-%H%M%S).tsv

{
    echo "#population GBE_ID n_finished_batch" | tr ' ' '\t'

    find /scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/*-batch/logs -type f -name "*log" \
    | grep -v 'QT_ALL' | sed -e 's%/%.%g' | awk -v FS='.' '{print $(NF-6), $(NF-2)}' | sed -e 's/-batch//g' | sort --parallel 2 | uniq -c \
    | sort --parallel 2 -k2,2 -k3,3V -k1,1 | awk -v OFS='\t' '{print $2, $3, $1}'

} > ${out_f}

echo ${out_f}

#!/bin/bash
set -beEuo pipefail

batch_idx=${SLURM_ARRAY_TASK_ID:=1}
jobs_prefix="7_ldsc_rg.20200706-144408"

if [ $# -gt 0 ] ; then batch_idx=$1 ; fi
if [ $# -gt 1 ] ; then jobs_prefix=$2 ; fi

jobs_ukb="${jobs_prefix}.ukb.tsv"
jobs_finngen="${jobs_prefix}.finngen.tsv"

n_ukb=$( cat ${jobs_ukb} | egrep -v '^#' | wc -l )

finngen_idx=$(perl -e "print(int((${batch_idx}-1)/${n_ukb})+1)")
ukb_idx=$(perl -e "print(${batch_idx}-((${finngen_idx}-1) * ${n_ukb}))")

FinnGen=$(   cat ${jobs_finngen} | egrep -v '^#' | awk -v nr=${finngen_idx} '(NR==nr){print $1}' )
GBE_ID=$(    cat ${jobs_ukb}     | egrep -v '^#' | awk -v nr=${ukb_idx}     '(NR==nr){print $1}' )
FinnGen_f=$( cat ${jobs_finngen} | egrep -v '^#' | awk -v nr=${finngen_idx} '(NR==nr){print $2}' )
WB_f=$(      cat ${jobs_ukb}     | egrep -v '^#' | awk -v nr=${ukb_idx}     '(NR==nr){print $2}' )

src="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ldsc_rg.sh"
finngen_d="/scratch/groups/mrivas/public_data/summary_stats/finngen_r3"

out_f="${finngen_d}/ldsc/UKB_WB_rg/${FinnGen}.${GBE_ID}"

if [ ! -d $(dirname ${out_f}) ] ; then mkdir -p $(dirname ${out_f}) ; fi

if [ ! -s ${out_f}.log ] ; then
    bash ${src} --scratch ${FinnGen_f} ${WB_f} ${out_f}
else
    echo ${out_f}.log
fi


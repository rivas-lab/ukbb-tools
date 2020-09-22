#!/bin/bash
set -beEuo pipefail

GBE_ID_1=$1
GBE_ID_2=$2
if [ $# -gt 2 ] ; then pop=$3 ; else pop="metal" ; fi
if [ $# -gt 3 ] ; then geno_dataset=$4 ; else geno_dataset="array-combined" ; fi

ldsc_dir="/oak/stanford/groups/mrivas/ukbb24983/${geno_dataset}/ldsc"
UKB_f1="${ldsc_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID_1}.${geno_dataset}.sumstats.gz"
UKB_f2="${ldsc_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID_2}.${geno_dataset}.sumstats.gz"
out_f="${ldsc_dir}/rg/${pop}.${GBE_ID_1}.${GBE_ID_2}"
src="/oak/stanford/groups/mrivas/users/${USER}/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ldsc_rg.sh"

if [ ! -d $(dirname ${out_f}) ] ; then mkdir -p $(dirname ${out_f}) ; fi

if [ ! -f ${out_f}.log ] ; then
    bash ${src} --scratch ${UKB_f1} ${UKB_f2} ${out_f}
fi


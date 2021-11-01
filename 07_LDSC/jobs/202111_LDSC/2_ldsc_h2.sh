#!/bin/bash
set -beEuo pipefail

GBE_ID=$1
if [ $# -gt 1 ] ; then pop=$2 ; else pop="white_british" ; fi
if [ $# -gt 2 ] ; then geno_dataset=$3 ; else geno_dataset="array-combined" ; fi

ldsc_dir="/oak/stanford/groups/mrivas/ukbb24983/${geno_dataset}/ldsc"
UKB_f="${ldsc_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}.${geno_dataset}.sumstats.gz"
out_f="${ldsc_dir}/h2/${pop}.${GBE_ID}"
src="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ldsc_h2.sh"

if [ ! -d $(dirname ${out_f}) ] ; then mkdir -p $(dirname ${out_f}) ; fi

if [ ! -f ${out_f}.log ] ; then
    bash ${src} --scratch ${UKB_f} ${out_f}
fi

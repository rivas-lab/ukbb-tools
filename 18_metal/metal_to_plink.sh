#!/bin/bash
set -beEuo pipefail

src="/oak/stanford/groups/mrivas/users/$USER/repos/rivas-lab/ukbb-tools/04_gwas/extras/sumstats_to_plink.sh"

# in_f="/oak/stanford/groups/mrivas/ukbb24983/array-combined/metal/20200717/HC382.metal.tsv.gz"
in_f=$1
out_f=$(dirname $(dirname ${in_f}))/$(basename $(dirname ${in_f}))_plink/$(basename ${in_f} .metal.tsv.gz).metal.plink.tsv.gz

if [ ! -d $(dirname ${out_f}) ] ; then mkdir -p $(dirname ${out_f}) ; fi

GBE_ID=$(basename ${in_f} .metal.tsv.gz)

GBE_CAT=$(echo ${GBE_ID} | sed -e "s/[0-9]//g")

if [ ! -f ${out_f} ] ; then

    if [ "${GBE_CAT}" == "QT_FC" ] || [ "${GBE_CAT}" == "INI" ] ; then
        bash ${src} ${in_f}
    else
        bash ${src} --logit ${in_f}
    fi | bgzip -l9 > ${out_f}

fi

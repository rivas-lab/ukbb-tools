#!/bin/bash
set -beEuo pipefail

GBE_ID=$1
pop=white_british

##################
GBE_ID_to_type () {
    # Input: GBE_ID
    # Output: linear or logistic
    local GBE_ID=$1
    # logistic: BIN|BIN_FC|cancer|FH|HC
    # linear: INI|QT_FC
    GBE_ID_prefix="$(echo ${GBE_ID} | sed -e 's/[0-9]*//g')"
    if [ "${GBE_ID_prefix}" == "INI" ] || [ "${GBE_ID_prefix}" == "QT_FC" ] ; then
        echo "linear"
    else
        echo "logistic"
    fi
}

##################
ukb_ldsc_munge_wrapper_sh="/oak/stanford/groups/mrivas/users/$USER/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ukb_ldsc_munge_wrapper.sh"
in_d=/oak/stanford/groups/mrivas/ukbb24983/cal/gwas/current/${pop}
out_d=/scratch/groups/mrivas/ukbb24983/array_combined/ldsc/current/${pop}
# # we have a copy in SCRATCH
# ld_scores=/scratch/groups/mrivas/projects/h2-estimation/private_output/ukb_ld_scores/TWINSUK
type=$(GBE_ID_to_type ${GBE_ID}) # "linear"
in_f=$(readlink -f ${in_d}/ukb24983_v2_hg19.${GBE_ID}.array-combined.glm.$(echo ${type} | sed -e "s/logistic/logistic.hybrid/g").gz)
out_f="${out_d}/ukb24983_v2_hg19.${GBE_ID}.array-combined.glm.sumstats.gz"

if [ ! -d "${out_d}" ] ; then mkdir -p ${out_d} ; fi

if [ ! -f "${out_f}" ] && [ ! -f "${out_f%.gz}.log" ] ; then
    echo bash ${ukb_ldsc_munge_wrapper_sh} ${out_f} ${GBE_ID} ${type} ${in_f}
    
    echo ""

    bash ${ukb_ldsc_munge_wrapper_sh} ${out_f} ${GBE_ID} ${type} ${in_f}

    #  \
    # ${ld_scores}
fi

#!/bin/bash
set -beEuo pipefail

GBE_ID=$1

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

src="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ldsc_munge.sh"
pop="white_british"
type=$(GBE_ID_to_type ${GBE_ID}) # "linear"
in_d=/oak/stanford/groups/mrivas/ukbb24983/cal/gwas/current/${pop}
out_d=/scratch/groups/mrivas/ukbb24983/array_combined/ldsc/202005/${pop}
in_f=$(readlink -f ${in_d}/ukb24983_v2_hg19.${GBE_ID}.array-combined.glm.$(echo ${type} | sed -e "s/logistic/logistic.hybrid/g").gz)
out_f="${out_d}/ukb24983_v2_hg19.${GBE_ID}.array-combined.glm.sumstats.gz"

if [ ! -d "${out_d}" ] ; then mkdir -p ${out_d} ; fi
if [ ! -f "${out_f}" ] && [ ! -f "${out_f%.sumstats.gz}.sumstats.gz" ] ; then bash ${src} --scratch ${in_f} ${out_f} ; fi


#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename ${SRCNAME})
VERSION="0.0.1"
NUM_POS_ARGS="1"

#####################

phe_to_gwas_location () {
    local phe_file=$1
    local geno_data=$2
    local pop=$3

    echo $(dirname $(dirname $phe_file))/${pop} | sed -e "s%phenotypedata%${geno_data}/gwas%g"
}

get_plink_suffix () {
    local GBE_ID=$1

    GBE_CAT=$(echo $GBE_ID | sed -e "s/[0-9]//g")

    if [ "${GBE_CAT}" == "QT_FC" ] || [ "${GBE_CAT}" == "INI" ] ; then
        echo glm.linear
    else
        echo glm.logistic.hybrid
    fi
}

#####################

pop=$1 # "others"

#####################

info_file=$(dirname ${SRCDIR})/05_gbe/phenotype_info.tsv
geno_data="array-combined"
symlink_dir=/oak/stanford/groups/mrivas/private_data/ukbb/24983/array-combined/gwas/current/${pop}
if [ ! -d ${symlink_dir} ] ; then mkdir -p ${symlink_dir} ; fi

cat ${info_file} | awk '{print $1, $NF}' | while read GBE_ID phe_file ; do
    new_basename="ukb24983_v2_hg19.${GBE_ID}.array-combined.$(get_plink_suffix $GBE_ID).gz"
    old_basename="ukb24983_v2_hg19.array-combined.${GBE_ID}.$(get_plink_suffix $GBE_ID).gz"

    if [ -f "${symlink_dir}/${old_basename}" ] ; then
        data_d=$(phe_to_gwas_location ${phe_file} ${geno_data} ${pop})
        if [ ! -d ${data_d} ] ; then mkdir -p ${data_d} ; fi
        mv ${data_d}/${old_basename} ${data_d}/${new_basename}
        rm ${symlink_dir}/${old_basename}
        ln -s ${data_d}/${new_basename} ${symlink_dir}
        echo ${symlink_dir}/${new_basename}
    fi
done

exit 0

# This is a script for a one-time fix of sumstats names
- wrong:   ukb24983_v2_hg19.array-combined.cancer1043.glm.logistic.hybrid.gz
- correct: ukb24983_v2_hg19.cancer1043.array-combined.glm.logistic.hybrid.gz

bash 20200606-gwas-additional-pops-mv.fix.sh others
bash 20200606-gwas-additional-pops-mv.fix.sh related

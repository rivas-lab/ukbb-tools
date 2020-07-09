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

get_sumstats_dir () {
    local pop=$1

    if [ "${pop}" == "others" ] ; then
        echo "/oak/stanford/groups/mrivas/projects/gwas_others"
    elif [ "${pop}" == "related" ] ; then
        echo "/oak/stanford/groups/mrivas/projects/related"
    else
        echo "unsupported pop ${pop}" >&2 ; exit 1
    fi
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
sumstats_dir=$(get_sumstats_dir $pop)
symlink_dir=/oak/stanford/groups/mrivas/private_data/ukbb/24983/array-combined/gwas/current/${pop}
if [ ! -d ${symlink_dir} ] ; then mkdir -p ${symlink_dir} ; fi

cat ${info_file} | awk '{print $1, $NF}' | while read GBE_ID phe_file ; do
    sumstats_f="${sumstats_dir}/ukb24983_v2_hg19.${GBE_ID}.array-combined.$(get_plink_suffix $GBE_ID).gz"
    if [ -f "${sumstats_f}" ] ; then
        data_d=$(phe_to_gwas_location ${phe_file} ${geno_data} ${pop})
        if [ ! -d ${data_d} ] ; then mkdir -p ${data_d} ; fi
        # mv file
        mv $sumstats_f $data_d/
        # create a symlink
        ln -s $data_d/$(basename $sumstats_f) $symlink_dir/
        echo $sumstats_f
    fi
done

exit 0
# this script move the sumstats in the array-combined/gwas directory based on the information in the 
# phenotype information file (05_gbe/phenotype_info.tsv) and generate symlink to the gwas/current dir
bash 20200606-gwas-additional-pops-mv.sh others
bash 20200606-gwas-additional-pops-mv.sh related

# /oak/stanford/groups/mrivas/private_data/ukbb/24983/array-combined/gwas/current/others
# /oak/stanford/groups/mrivas/private_data/ukbb/24983/array-combined/gwas/current/related

# 2020/6/16 11:50 --> YT applied this script. We need to run it again once we have the results from the second half of the GWAS
#  So the ToDo items are:
#  - wait and let it finish the
for pop in others related ; do bash 20200606-gwas-additional-pops-mv.sh $pop ; done | tee 20200606-gwas-additional-pops-mv.20200626.log

for pop in others related ; do bash 20200606-gwas-additional-pops-mv.sh $pop ; done | tee 20200606-gwas-additional-pops-mv.20200627.log

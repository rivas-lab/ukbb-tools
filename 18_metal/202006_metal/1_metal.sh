#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename ${SRCNAME})
VERSION="0.0.1"
NUM_POS_ARGS="1"

#####################

get_plink_suffix () {
    local GBE_ID=$1

    GBE_CAT=$(echo $GBE_ID | sed -e "s/[0-9]//g")

    if [ "${GBE_CAT}" == "QT_FC" ] || [ "${GBE_CAT}" == "INI" ] ; then
        echo glm.linear
    else
        echo glm.logistic.hybrid
    fi
}

############################################################
# tmp dir
############################################################
tmp_dir_root="$LOCAL_SCRATCH"
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT

#####################

GBE_ID=$1
# GBE_ID=INI50

#####################

info_file=$(dirname $(dirname ${SRCDIR}))/05_gbe/phenotype_info.tsv
metal_script=$(dirname ${SRCDIR})/run_metal.sh
geno_data="array-combined"
pops=(white_british non_british_white african s_asian e_asian related others)
data_dir=/oak/stanford/groups/mrivas/private_data/ukbb/24983/${geno_data}
metal_dir=${data_dir}/metal/20200616
in_file_list=${tmp_dir}/metal.input.lst

#####################

if [ ! -d ${metal_dir} ] ; then mkdir -p ${metal_dir} ; fi

cat ${info_file} | awk -v GBE_ID=${GBE_ID} '($1 == GBE_ID){print $1, $NF}' | while read GBE_ID phe_file ; do
    for pop in ${pops[@]} ; do
        sumstats_l=$data_dir/gwas/current/${pop}/ukb24983_v2_hg19.${GBE_ID}.${geno_data}.$(get_plink_suffix ${GBE_ID}).gz
        if [ -f "${sumstats_l}" ] ; then
            readlink -f ${sumstats_l}
        fi
    done
done > ${in_file_list}

n_in_files=$(cat ${in_file_list} | wc -l)

if [ "${n_in_files}" -gt 0 ] ; then
    bash ${metal_script} -o ${metal_dir}/${GBE_ID} -f ${in_file_list}
fi

exit 0
#######################
usage:
ml load R/3.6 gcc/6 # snpnet_yt # or your favorite R env
bash 1_metal.sh INI50
bash 1_metal.sh HC276

#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename ${SRCNAME})
VERSION="0.1.0"
NUM_POS_ARGS="1"

############################################################
# constants
############################################################

info_file=$(dirname $(dirname ${SRCDIR}))/05_gbe/phenotype_info.tsv
metal_script=$(dirname ${SRCDIR})/run_metal.sh

############################################################
# functions
############################################################

show_default_helper () {
    cat ${SRCNAME} | grep -n Default | tail -n+3 | awk -v FS=':' '{print $1}' | tr "\n" "\t"
}

show_default () {
    cat ${SRCNAME} \
        | tail -n+$(show_default_helper | awk -v FS='\t' '{print $1+1}') \
        | head  -n$(show_default_helper | awk -v FS='\t' '{print $2-$1-1}')
}

usage () {
cat <<- EOF
	$PROGNAME (version $VERSION)
	Run Metal for the UKB sumstats

	Usage: $PROGNAME [options] GBE_ID
	  GBE_ID        The GBE_ID

	Options:

	Default configurations:
	  info_file=${info_file}
	  metal_script=${metal_script}
EOF
    show_default | awk -v spacer="  " '{print spacer $0}'
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

############################################################
# tmp dir
############################################################
tmp_dir_root="$LOCAL_SCRATCH"
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
# echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT

############################################################
# parser start
############################################################
## == Default parameters (start) == ##
ukb_data_dir=/oak/stanford/groups/mrivas/private_data/ukbb/24983
geno_data=array-combined
pops_str=white_british,non_british_white,african,s_asian,e_asian,related,others
metal_freeze_v=20200717
## == Default parameters (end) == ##

declare -a params=()
for OPT in "$@" ; do
    case "$OPT" in
        '-h' | '--help' )
            usage >&2 ; exit 0 ;
            ;;
        '-v' | '--version' )
            echo $VERSION ; exit 0 ;
            ;;
        '--ukb_data_dir' )
            ukb_data_dir=$2 ; shift 2 ;
            ;;
        '--geno_data' )
            geno_data=$2 ; shift 2 ;
            ;;
        '--pops_str' )
            pops_str=$2 ; shift 2 ;
            ;;
        '--metal_freeze_v' )
            metal_freeze_v=$2 ; shift 2 ;
            ;;
        '--'|'-' )
            shift 1 ; params+=( "$@" ) ; break
            ;;
        -*)
            echo "$PROGNAME: illegal option -- '$(echo $1 | sed 's/^-*//')'" 1>&2 ; exit 1
            ;;
        *)
            if [[ $# -gt 0 ]] && [[ ! "$1" =~ ^-+ ]]; then
                params+=( "$1" ) ; shift 1
            fi
            ;;
    esac
done

if [ ${#params[@]} -lt ${NUM_POS_ARGS} ]; then
    echo "${PROGNAME}: ${NUM_POS_ARGS} positional arguments are required" >&2
    usage >&2 ; exit 1 ;
fi

GBE_ID="${params[0]}"

############################################################

# check R environmnet
Rscript /dev/stdin << EOF
suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))
EOF

metal_dir=${ukb_data_dir}/${geno_data}/metal/${metal_freeze_v}
in_file_list=${tmp_dir}/metal.input.lst
if [ ! -d ${metal_dir} ] ; then mkdir -p ${metal_dir} ; fi

echo ${pops_str} | tr ',' '\n' | while read pop ; do
    sumstats_l=${ukb_data_dir}/${geno_data}/gwas/current/${pop}/ukb24983_v2_hg19.${GBE_ID}.${geno_data}.$(get_plink_suffix ${GBE_ID}).gz
    if [ -f "${sumstats_l}" ] ; then
        readlink -f ${sumstats_l}
    else
        echo ${sumstats_l} >&2
    fi
done > ${in_file_list}

n_in_files=$(cat ${in_file_list} | wc -l)

if [ "${n_in_files}" -gt 0 ] ; then
    bash ${metal_script} -o ${metal_dir}/${GBE_ID} -f ${in_file_list}
else
    echo "We don't have sumstats!"
fi

exit 0
#######################
usage:
ml load R/3.6 gcc/6 # snpnet_yt # or your favorite R env
bash 1_metal.sh INI50
bash 1_metal.sh HC276

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
ukb_data_dir=/oak/stanford/groups/mrivas/private_data/ukbb/24983
geno_data=array-combined
metal_freeze_v=20200717
N_thr=100

############################################################
# body
############################################################

cat ${info_file} | egrep -v '^#' | awk -v FS='\t' -v thr=$N_thr '($7 >= thr && length($1)>0){print $1}' | while read GBE_ID ; do
    metal_prefix=${ukb_data_dir}/${geno_data}/metal/${metal_freeze_v}/${GBE_ID}.metal
    if [ ! -f ${metal_prefix}.tsv.gz ] && [ ! -f ${metal_prefix}.info.txt ] ; then
        echo ${GBE_ID}
    fi
done

#!/bin/bash
set -beEuo pipefail

LDSC_rg_log_dir="/oak/stanford/groups/mrivas/ukbb24983/array-combined/ldsc/rg/"

pop='metal' # population 
if [ $# -gt 0 ] ; then pop=$1 ; fi

LDSC_rg_f=$(dirname ${LDSC_rg_log_dir})/rg.${pop}.$(date +%Y%m%d-%H%M%S).tsv
if [ $# -gt 1 ] ; then LDSC_rg_f=$2 ; fi

############################################################
# tmp dir
############################################################
tmp_dir_root="${LOCAL_SCRATCH:=/tmp}"
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
# echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT
############################################################

############################################################
# 1) list the LDSC rg log files
############################################################
LDSC_rg_logs_f=${tmp_dir}/ldsc_rg.lst
find ${LDSC_rg_log_dir} -type f -name "${pop}.*.log" | sort -V > ${LDSC_rg_logs_f}

############################################################
# 2) aggregate the results
############################################################

echo "We found $(cat ${LDSC_rg_logs_f} | wc -l) files. Writing the results table to ${LDSC_rg_f}"

bash /oak/stanford/groups/mrivas/users/${USER}/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ldsc_rg_view.sh -l ${LDSC_rg_logs_f} \
| sed -e "s%/scratch%/oak/stanford%g" \
| sed -e "s%/oak/stanford/groups/mrivas/ukbb24983/array-combined/ldsc/metal/ukb24983_v2_hg19.%%g" \
| sed -e "s%.array-combined.sumstats.gz%%g" \
| awk -v FS='\t' -v pop=${pop} '{print pop, $0}' \
| sed -e "s/^p1/#population p1/g" | tr ' ' '\t'  > ${LDSC_rg_f}


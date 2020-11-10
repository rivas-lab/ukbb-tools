#!/bin/bash
set -beEuo pipefail

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

ml load R/3.6 gcc

data_d=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l
tmp_f=${tmp_dir}/wc_l.tsv
out_f=$(dirname ${data_d})/wc_l.$(date +%Y%m%d-%H%M%S).tsv

{
echo '#file NF_wc_l wc_l non_NA_wc_l' | tr ' ' '\t'

find ${data_d} -type f | sort -V | parallel -k 'cat {} | egrep -v "^#"'
} > ${tmp_f}

Rscript 7c_wc_l_aggregate.sub.R ${tmp_f} ${out_f} 

echo ${out_f}

#!/bin/bash
set -beEuo pipefail

GBE_ID=$1

jn=gwas_GBE_ID_${GBE_ID}

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

sbatch --mail-user=ytanigaw@stanford.edu --mail-type=END,FAIL \
--account=mrivas -p nih_s10 \
--nodes=1 --mem=13000 --cores=4 --time=2:00:00 \
--job-name=${jn} --output=logs/${jn}.%A_%a.out --error=logs/${jn}.%A_%a.err \
--array=1-100 \
1_plink.gwas.GBE_ID.scg.sh ${GBE_ID} | tee ${tmp_dir}/job.log

scontrol update job $(cat ${tmp_dir}/job.log | awk '{print $NF}') partition=nih_s10,batch 


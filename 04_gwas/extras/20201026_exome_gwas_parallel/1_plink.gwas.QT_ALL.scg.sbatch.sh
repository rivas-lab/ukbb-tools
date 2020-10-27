#!/bin/bash
set -beEuo pipefail

pop="white_british"

if [ $# -gt 0 ] ; then pop=$1 ; fi

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

# jn=gwas_QT_WB
jn=gwas_QT_${pop}

sbatch --mail-user=ytanigaw@stanford.edu --mail-type=END,FAIL \
--account=mrivas -p nih_s10 \
--nodes=1 --mem=40000 --cores=8 \
--time=2:00:00 \
--array=1-100 \
--job-name=${jn} --output=logs/${jn}.%A_%a.out --error=logs/${jn}.%A_%a.err \
1_plink.gwas.QT_ALL.scg.sh ${pop} | tee ${tmp_dir}/job.log

scontrol update job $(cat ${tmp_dir}/job.log | awk '{print $NF}') partition=nih_s10,batch

exit 0
# for WB
--time=12:00:00 \
--array=251-1000
--array=6,38,40,41,42,43,44,45,48,49,52,53,54,55,78,79,80,82,86,87,89,95,96,97,127,129,130,131,132,133,134,136,146,147,149,151,152,156,157


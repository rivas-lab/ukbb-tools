#!/bin/bash
set -beEuo pipefail

list=$1
if [ $# -gt 1 ] ; then offset=$2 ; else offset=0 ; fi

log_d=logs_20201006

if [ ! -d ${log_d} ] ; then mkdir -p ${log_d} ; fi

ml load resbatch

# sbatch -p mrivas,normal,owners \
sbatch -p mrivas --qos=high_p \
--time=1:0:00 --mem=8000 --nodes=1 --cores=1 --job-name=vep --output=${log_d}/vep.%A_%a.out --error=${log_d}/vep.%A_%a.err \
--array=1-933bank \
${parallel_sbatch_no_err_check_sh} \
4-4_resubmit.sh ${parallel_idx} 1 ${list} ${offset}

exit 0
/oak/stanford/groups/mrivas/ukbb24983/cnv/pgen/cnv.pvar
# 275181 lines, 275180 variants
# 98504 lines, 98503 variants
# split into 9851 pieces with up to 10 lines each
# 982 array jobs with 10 pieces each


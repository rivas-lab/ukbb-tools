#!/bin/bash
set -beEuo pipefail

GBE_lst=$1
if [ $# -gt 1 ] ; then batch_idx_s=$2 ; else batch_idx_s=11 ; fi

jn=gwas.WB.$(basename ${GBE_lst%.lst})

sbatch -p mrivas,normal,owners \
    --time=23:55:00 --mem=9000 --nodes=1 --cores=2 \
    --job-name=$jn --output=logs_scratch/$jn.%A_%a.out --error=logs_scratch/$jn.%A_%a.err \
    --array=1-$(cat $GBE_lst | egrep -v '#' | wc -l) ${parallel_sbatch_no_err_check_sh} \
    1_plink.gwas.GBE_ID.sh2-WB.oak.bugfixed.job.sh \
    ${GBE_lst} 1 ${batch_idx_s}


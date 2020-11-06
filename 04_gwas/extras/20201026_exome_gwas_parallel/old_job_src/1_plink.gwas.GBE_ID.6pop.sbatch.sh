#!/bin/bash
set -beEuo pipefail

GBE_lst=$1

jn=gwas.$(basename ${GBE_lst%.lst})

sbatch -p mrivas,normal,owners \
    --time=1-0:0:00 --mem=9000 --nodes=1 --cores=2 \
    --job-name=$jn --output=logs_scratch/$jn.%A_%a.out --error=logs_scratch/$jn.%A_%a.err \
    --array=1-$(cat $GBE_lst | egrep -v '#' | wc -l) ${parallel_sbatch_no_err_check_sh} \
    1_plink.gwas.GBE_ID.6pop.job.sh \
    ${GBE_lst} 1


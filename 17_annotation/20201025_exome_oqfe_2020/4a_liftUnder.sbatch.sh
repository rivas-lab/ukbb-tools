#!/bin/bash
set -beEuo pipefail

log_d=logs

if [ ! -d ${log_d} ] ; then mkdir -p ${log_d} ; fi

ml load resbatch

sbatch -p mrivas --qos=high_p --time=1:0:00 --mem=6000 --nodes=1 --cores=1 \
    --job-name=liftUnder --output=${log_d}/liftUnder.%A_%a.out --error=${log_d}/liftUnder.%A_%a.err \
    --array=2-900 ${parallel_sbatch_sh} 4a_liftUnder.sh ${parallel_idx} 1


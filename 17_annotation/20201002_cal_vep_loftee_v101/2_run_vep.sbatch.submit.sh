#!/bin/bash
set -beEuo pipefail

log_d=logs

if [ ! -d ${log_d} ] ; then mkdir -p ${log_d} ; fi

ml load resbatch

sbatch -p mrivas --qos=high_p --time=1:0:00 --mem=6000 --nodes=1 --cores=1 --job-name=vep --output=${log_d}/vep.%A_%a.out --error=${log_d}/vep.%A_%a.err --array=1-403 ${parallel_sbatch_sh} 2_run_vep.sh ${parallel_idx} 1
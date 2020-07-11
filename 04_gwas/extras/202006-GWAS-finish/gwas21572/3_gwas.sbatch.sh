#!/bin/bash
set -beEuo pipefail

batch_idx=${SLURM_ARRAY_TASK_ID:=1}
if [ $# -gt 0 ] ; then batch_idx=$1 ; fi

jobs="job.20200710-235051.tsv"

GBE_ID=$( cat ${jobs} | egrep -v '^#' | awk -v nr=${batch_idx} '(NR==nr){print $1}' )
pop=$(    cat ${jobs} | egrep -v '^#' | awk -v nr=${batch_idx} '(NR==nr){print $2}' )

echo bash 3_gwas.sh $GBE_ID $pop
bash 3_gwas.sh $GBE_ID $pop

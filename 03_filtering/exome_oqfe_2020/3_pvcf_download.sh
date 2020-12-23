#!/bin/bash
set -beEuo pipefail

if [ $# -gt 0 ] ; then
    idx=$1
else
    idx=${SLURM_ARRAY_TASK_ID:=1}
fi

echo ${idx}

data_d="/scratch/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/download"
idx_f=${data_d}/pvcf_blocks.txt

if [ "${idx}" -gt 0 ] && [ "${idx}" -lt 978 ]  && [ "${idx}" -ne 977 ] && [ "${idx}" -ne 953 ] ; then

    c=$(cat ${idx_f} | awk -v idx=${idx} '($1 == idx){print $2}')
    b=$(cat ${idx_f} | awk -v idx=${idx} '($1 == idx){print $3}')

    cd ${data_d}

    echo "./gfetch 23156 -c${c} -b${b}"
    ./gfetch 23156 -c${c} -b${b}

fi

exit 0
#####################

We tested the gfetch commands

3_pvcf_download.sh 953
3_pvcf_download.sh 977

We then submit the scripts to SLURM.

sbatch -p mrivas,normal --time=4:0:00 --mem=8000 --nodes=1 --cores=1 --job-name=pvcf_DL --output=logs/pvcf_DL.%A_%a.out --error=logs/pvcf_DL.%A_%a.err --array=1-976%5 ${parallel_sbatch_sh} 3_pvcf_download.sh ${parallel_idx} 1


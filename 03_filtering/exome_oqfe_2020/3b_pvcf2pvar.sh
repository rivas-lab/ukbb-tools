#!/bin/bash
set -beEuo pipefail

if [ $# -gt 0 ] ; then
    idx=$1
else
    idx=${SLURM_ARRAY_TASK_ID:=1}
fi

ml load plink2

echo ${idx}

data_d="/scratch/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/download"
idx_f=${data_d}/pvcf_blocks.txt

if [ "${idx}" -gt 0 ] && [ "${idx}" -lt 978 ] ; then

    c=$(cat ${idx_f} | awk -v idx=${idx} '($1 == idx){print $2}' | sed -e 's/23/X/g' | sed -e 's/24/Y/g')
    b=$(cat ${idx_f} | awk -v idx=${idx} '($1 == idx){print $3}')
    out="$(dirname ${data_d})/pvcf_pvar/ukb23156_c${c}_b${b}_v1"

    if [ ! -s ${out}.pvar.log ] ; then
        plink2 --threads 1 --memory 8000 --make-just-pvar zs --out ${out} \
            --vcf ${data_d}/$(basename ${out}).vcf.gz

        mv ${out}.log ${out}.pvar.log
    fi
fi

exit 0
#####################

We then submit the scripts to SLURM.

sbatch -p mrivas,normal --time=30:00 --mem=8000 --nodes=1 --cores=1 --job-name=pvcf2pvar --output=logs/pvcf2pvar.%A_%a.out --error=logs/pvcf2pvar.%A_%a.err --array=2-976 ${parallel_sbatch_sh} 3b_pvcf2pvar.sh ${parallel_idx} 1


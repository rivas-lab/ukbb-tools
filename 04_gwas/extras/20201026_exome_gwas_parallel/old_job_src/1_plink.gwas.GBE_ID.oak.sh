#!/bin/bash
set -beEuo pipefail

batch_idx=${SLURM_ARRAY_TASK_ID:=1}
pop=$1
GBE_ID=$2

if [ $# -gt 2 ] ; then batch_idx=$3 ; fi

out_d=/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_scg/${pop}-batch
#out_d=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/${pop}-batch

bash 1_plink.gwas.sh \
--GBE_ID ${GBE_ID} \
--pop ${pop} \
--n_batch $([ "${pop}" == "white_british" ] && echo "100" || echo "100") \
--pfile /scratch/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/ukb24983_exomeOQFE \
${batch_idx} \
${out_d}

exit 0
--pfile /oak/stanford/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/ukb24983_exomeOQFE \

sbatch -p mrivas --qos=high_p --nodes=2 --mem=9000 --cores=1 \
    --time=2:00:00 --job-name=gwas \
    --output=logs_scratch/gwas.%A_%a.out --error=logs_scratch/gwas.%A_%a.err \
    1_plink.gwas.GBE_ID.oak.sh white_british HC269


#!/bin/bash
set -beEuo pipefail

batch_idx=${SLURM_ARRAY_TASK_ID:=1}
pop=$1

if [ $# -gt 1 ] ; then batch_idx=$2 ; fi

#out_d=/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/${pop}-batch
#out_d=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/${pop}-batch
out_d=/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_scg/${pop}-batch

ml load zstd

bash 1_plink.gwas.sh \
--cores 8 \
--mem 32000 \
--pop ${pop} \
--n_batch $([ "${pop}" == "white_british" ] && echo "1000" || echo "100") \
--pfile /oak/stanford/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/ukb24983_exomeOQFE \
--QT_ALL \
${batch_idx} \
${out_d}


#!/bin/bash
set -beEuo pipefail

batch_idx=${SLURM_ARRAY_TASK_ID:=1}
GBE_ID=$1

if [ $# -gt 1 ] ; then batch_idx=$2 ; fi

pop=white_british
out_d=/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_scg/${pop}-batch

ml load zstd

bash 1_plink.gwas.sh \
--cores 4 \
--mem 16000 \
--GBE_ID ${GBE_ID} \
--pop ${pop} \
--pfile /oak/stanford/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/ukb24983_exomeOQFE \
${batch_idx} \
${out_d}


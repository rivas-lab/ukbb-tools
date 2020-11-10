#!/bin/bash
set -beEuo pipefail

batch_idxs=${SLURM_ARRAY_TASK_ID:=1}
pop=$1
GBE_ID=$2
if [ $# -gt 2 ] && [ $3 == "oak" ] ; then oak="TRUE" ; else oak="FALSE" ; fi
if [ $# -gt 3 ] ; then batch_idxs=$4 ; fi

if [ "${oak}" == "TRUE" ] ; then
    out_d=/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_scg/${pop}-batch
else
    out_d=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/${pop}-batch
fi

for batch_idx in $(echo ${batch_idxs} | tr ',' ' ') ; do

echo $batch_idx

! bash 1_plink.gwas.sh \
--GBE_ID ${GBE_ID} \
--pop ${pop} \
--pfile /scratch/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/ukb24983_exomeOQFE \
${batch_idx} \
${out_d}

done

exit 0
###############################################

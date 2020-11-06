#!/bin/bash
set -beEuo pipefail

batch_idx=${SLURM_ARRAY_TASK_ID:=1}
pop=$1
GBE_ID=$2

if [ $# -gt 2 ] ; then batch_idx=$3 ; fi

#out_d=/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/${pop}-batch
out_d=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/${pop}-batch

bash 1_plink.gwas.sh \
--GBE_ID ${GBE_ID} \
--pop ${pop} \
--n_batch $([ "${pop}" == "white_british" ] && echo "100" || echo "100") \
--pfile /scratch/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/ukb24983_exomeOQFE \
${batch_idx} \
${out_d}

exit 0


--pfile /oak/stanford/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/ukb24983_exomeOQFE \
Start time: Sun Oct 25 23:03:41 2020
End time: Sun Oct 25 23:11:17 2020
QT_ALL 8 min for 1/1000 batch (10449 variants) for 50k (~34k)
HC269 less than 1 min for 1/1000 batch (10449 variants) for 50k (~34k)
HC382 less than 1 min for 1/1000 batch (10449 variants) for 50k (~34k)
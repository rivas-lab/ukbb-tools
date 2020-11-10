#!/bin/bash
set -beEuo pipefail

batch_idx=${SLURM_ARRAY_TASK_ID:=1}

out_d=/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002-dev

bash 1_plink.gwas.sh \
--pfile /oak/stanford/groups/mrivas/ukbb24983/exome/pgen/ukb24983_exome \
--n_batch 1000 \
--overwrite \
--QT_ALL \
--GBE_ID HC269 \
${batch_idx} \
${out_d}

exit 0
Start time: Sun Oct 25 23:03:41 2020
End time: Sun Oct 25 23:11:17 2020
QT_ALL 8 min for 1/1000 batch (10449 variants) for 50k (~34k)
HC269 less than 1 min for 1/1000 batch (10449 variants) for 50k (~34k)
HC382 less than 1 min for 1/1000 batch (10449 variants) for 50k (~34k)

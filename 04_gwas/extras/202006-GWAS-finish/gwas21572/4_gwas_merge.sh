#!/bin/bash
set -beEuo pipefail

batch_idx=${SLURM_ARRAY_TASK_ID:=1}
if [ $# -gt 0 ] ; then batch_idx=$1 ; fi

ml load R/3.6 gcc

jobs="4_combined.idx.tsv"

symlink=$( cat ${jobs} | egrep -v '^#' | awk -v nr=${batch_idx} '(NR==nr){print $3}' )
add_res=$( cat ${jobs} | egrep -v '^#' | awk -v nr=${batch_idx} '(NR==nr){print $4}' )

bash /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/extras/202006-GWAS-finish/gwas369/4_combine.sh $symlink $add_res

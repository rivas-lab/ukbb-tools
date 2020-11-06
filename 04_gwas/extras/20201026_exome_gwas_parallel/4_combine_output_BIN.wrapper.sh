#!/bin/bash
set -beEuo pipefail

job_idx=$1
index_file=$2

pop=$(    cat ${index_file} | egrep -v '^#' | awk -v nr=${job_idx} '(NR == nr){print $1}' )
GBE_ID=$( cat ${index_file} | egrep -v '^#' | awk -v nr=${job_idx} '(NR == nr){print $2}' )

bash 4_combine_output_BIN.sh ${pop} ${GBE_ID}

exit 0
##########################
date ; bash 4_combine_output_BIN.wrapper.sh 1 3b_merge_job_list.20201028-194315.tsv ; date
2 min

sbatch -p mrivas,normal,owners --time=1:0:00 --mem=6000 --nodes=1 --cores=1 \
--job-name=combine --output=logs_scratch/combine.%A_%a.out --error=logs_scratch/combine.%A_%a.err \
--array=1-753 ${parallel_sbatch_sh} 4_combine_output_BIN.wrapper.sh ${parallel_idx} 10 3b_merge_job_list.20201028-214235.tsv

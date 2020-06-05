#!/bin/bash

ml load snpnet_yt # this is Yosuke's R env
ml load htslib

batch_idx=${SLURM_ARRAY_TASK_ID:=724}
batch_size=20
batch_job_file="202005_combine_chrX_check.20200604-220539.lst"
start_idx=$( perl -e "print(int(((${batch_idx} - 1) * ${batch_size}) + 1))" )
end_idx=$(   perl -e "print(int(  ${batch_idx}      * ${batch_size}))" )

cat ${batch_job_file} | awk -v start_idx=${start_idx} -v end_idx=${end_idx} 'start_idx <= NR && NR <= end_idx' |
while read pop GBE_ID ; do

    echo $pop $GBE_ID
    bash 202005_combine_chrX.sh $pop $GBE_ID

done


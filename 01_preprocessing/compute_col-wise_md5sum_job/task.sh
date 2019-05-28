#!/bin/bash
set -beEu -o pipefail

PROGNAME=$(basename $(readlink -f $0))
VERSION="1.0"

# pass the following paramters to your script
task_id=$1

job_tsv_file="jobs.tsv"
batch_size=50
repo_dir="/oak/stanford/groups/mrivas/users/$USER/repos/rivas-lab/ukbb-tools/"
script="${repo_dir}/01_preprocessing/compute_col-wise_md5sum_ary_job_task.sh"
out_dir="${repo_dir}/01_preprocessing/compute_col-wise_md5sum_job/out"

get_line_from_job_tsv () {
	job_tsv=$1
	task_id=$2
	
	cat $job_tsv | egrep -v '^#' | awk -v row=$task_id 'NR == row'
}

# create a temp directory
tmp_root=$LOCAL_SCRATCH
tmp_dir="$(mktemp -d -p ${tmp_root} tmp-${PROGNAME}-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
handler_exit () { rm -rf ${tmp_dir} ; }
trap handler_exit EXIT

out_file=${out_dir}/${task_id}.tsv
tmp_out_file=${tmp_dir}/$(basename ${out_file})

tab_file=$(get_line_from_job_tsv ${job_tsv_file} ${task_id} | awk '{print $2}')
per_tab_job_idx=$(get_line_from_job_tsv ${job_tsv_file} ${task_id} | awk '{print $1}')

if [ $task_id -eq 1 ] ; then
    bash ${script} ${tab_file} ${per_tab_job_idx} ${batch_size} -H > ${tmp_out_file}
else
    bash ${script} ${tab_file} ${per_tab_job_idx} ${batch_size} > ${tmp_out_file}
fi

cp ${tmp_out_file} ${out_file}


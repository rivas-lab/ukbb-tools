#!/bin/bash
set -beEu -o pipefail

# pass the following paramters to your script
if [ $# -gt 0 ] ; then job_task_id=${1} ; else job_task_id=1 ; fi

my_repos_dir="$OAK/users/$USER/repos"
repo_dir="${my_repos_dir}/rivas-lab/ukbb-tools"
script="${repo_dir}/07_LDSC/scripts/ukb_ldsc_munge_rg.sh"

job_tsv=$(dirname $0)/jobs.tsv
dataset_name="scz_autism"

get_line_from_job_tsv () {
    job_tsv=$1
    task_id=$2	
    cat $job_tsv | egrep -v '^#' | awk -v row=$task_id 'NR == row'
}

run_task () {	
    task_id=$1
    sumstats1="$(get_line_from_job_tsv $job_tsv ${task_id} | awk -v FS='\t' '{print $1}')"
    sumstats2="$(get_line_from_job_tsv $job_tsv ${task_id} | awk -v FS='\t' '{print $2}')"
    echo bash ${script} ${sumstats1} ${sumstats2} ${dataset_name}
    if [ ! -z "${sumstats1}" ] && [ ! -z "${sumstats2}" ] && [ ${sumstats1} != "" ] && [ ${sumstats2} != "" ] ; then
        bash ${script} ${sumstats1} ${sumstats2} ${dataset_name}
    fi
}

# For each task in the array job, we will run the script 10 times with different parameters.
# For example, $job_task_id == 1 corresponds to task_id = 1, 2, ... , 10

for offset in $(seq 1 10) ; do
    task_id=$(perl -e "print(10 * ${job_task_id} + ${offset} - 10)")
	#task_id=$(expr "${job_task_id}0" - $offset + 1)
	echo "----------" >&2
	echo "processing $task_id" >&2
	run_task $task_id
	echo "" >&2
done


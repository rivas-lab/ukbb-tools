#!/bin/bash
set -beEu -o pipefail

# pass the following paramters to your script
job_task_id=$1

my_repos_dir="$OAK/users/$USER/repos"
repo_dir="${my_repos_dir}/rivas-lab/ukbb-tools"
script="${repo_dir}/07_LDSC/scripts/ukb_ldsc_munge_h2.sh"

run_task () {
	task_id=$1
	
	sumstats=$(cat $(dirname $0)/jobs.tsv | awk -v row=$task_id '(NR == row){print $1}')

	if [ $(echo $sumstats | wc -c) -lt 2 ] ; then 
		echo "specified task_id ($task_id) is out of the expected range" >&2 
	else
		bash $script ${sumstats}
	fi

}

# For each task in the array job, we will run the script 10 times with different parameters.
# For example, $job_task_id == 1 corresponds to task_id = 1, 2, ... , 10

for offset in $(seq 1 10) ; do
	task_id=$(expr "${job_task_id}0" - $offset + 1)
	echo "----------" >&2
	echo "processing $task_id" >&2
	run_task $task_id
	echo "" >&2
done


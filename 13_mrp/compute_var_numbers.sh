#!/bin/bash
#SBATCH --job-name=make_R_phen_row
#SBATCH --output=mrp_logs/make_R_phen_row.%A_%a.out
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=16000
#SBATCH --time=08:00:00
#SBATCH -p normal,owners,mrivas

# define functions
usage () {
    echo "$0: Script that estimates one row of the R_phen matrix."
    echo "usage: sbatch --array=1-<number of rows to calculate> $0 start_idx (inclusive)"
    echo "e.g. sbatch --array=1-1000 $0 1"
}

# get core and memory settings from the header -- passed to gwas script below
cores=$(              cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(                cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )
log_dir=$( dirname $( cat $0 | egrep '^#SBATCH --output=' | awk -v FS='=' '{print $NF}' ))

# check number of command line args and dump usage if that's not right
if [ $# -ne 1 ] ; then usage >&1 ; exit 1 ; fi

# load sofware, dump which versions are used
export MODULEPATH="/home/groups/mrivas/.modules:${MODULEPATH}"
ml load python/2.7.13

# job start header (for use with array-job module)
_SLURM_JOBID=${SLURM_JOBID:=0} # use 0 for default value (for debugging purpose)
_SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:=1}
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname = $(hostname) SLURM_JOBID = ${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${_SLURM_ARRAY_TASK_ID}" >&1

# get phenotypes to run
start_idx=$1
this_idx=$_SLURM_ARRAY_TASK_ID
cur_idx=` echo "$start_idx + $this_idx - 1" | bc `

GBE_ID=$(cat sumstat_paths.tsv | grep -v path | awk -v start_idx=$start_idx -v this_idx=$this_idx 'NR==(start_idx + this_idx - 1) {print $1}' )

echo $GBE_ID
echo $cur_idx

/share/software/user/open/python/2.7.13/bin/python compute_var_numbers.py $GBE_ID $cur_idx

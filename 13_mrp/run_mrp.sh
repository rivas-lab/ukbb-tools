#!/bin/bash
#SBATCH --job-name=MRP
#SBATCH --output=mrp_logs/run_mrp.%A_%a.out
#SBATCH  --error=mrp_logs/run_mrp.%A_%a.err
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=25600
#SBATCH --time=1-00:00:00
#SBATCH -p normal,owners,mrivas

set -beEuo pipefail

# define functions
usage () {
    echo "$0: MRP re-run script for the exome data"
    echo "usage: sbatch --array=1-<number of array jobs> $0 start_idx (inclusive.)"
    echo ''
    echo '  You may check the status of the job (which jobs are finished) using the array-job module:'
    echo '  $ ml load array-job'
    echo '  $ array-job-find_ok.sh rerun_logs'
}

software_versions () {
    which plink2
    which bgzip
}

# get core and memory settings from the header -- passed to gwas script below
cores=$(              cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(                cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )
log_dir=$( dirname $( cat $0 | egrep '^#SBATCH --output=' | awk -v FS='=' '{print $NF}' ))

# check number of command line args and dump usage if that's not right
if [ $# -lt 1 ] ; then usage >&2 ; exit 1 ; fi

# load sofware, dump which versions are used
export MODULEPATH="/home/groups/mrivas/.modules:${MODULEPATH}"
ml load htslib
ml load python/3.6.1
if grep -q "CPU_GEN:HSW\|CPU_GEN:BDW\|CPU_GEN:SKX" <(a=$(hostname); sinfo -N -n ${a::-4} --format "%50f"); then
   # AVX2 is suitable for use on this node if CPU is recent enough
   ml load plink2/20190402
else
   ml load plink2/20190402-non-AVX2
fi

software_versions >&2

# job start header (for use with array-job module)
_SLURM_JOBID=${SLURM_JOBID:=0} # use 0 for default value (for debugging purpose)
_SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:=1}
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname = $(hostname) SLURM_JOBID = ${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${_SLURM_ARRAY_TASK_ID}" >&2

# get phenotypes to run
start_idx=$1
this_idx=$_SLURM_ARRAY_TASK_ID

min_N_count=10
GBE_ID=$(cat ../05_gbe/phenotype_info.tsv | awk -v min_N=${min_N_count} 'NR > 1 && $7 >= min_N' | egrep -v MED | awk -v start_idx=$start_idx -v this_idx=$this_idx 'NR==(start_idx + this_idx - 1) {print $1}' )

python3 mrp.py $GBE_ID

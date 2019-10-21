#!/bin/bash
#SBATCH --job-name=bgen2plink
#SBATCH --output=run_logs/run_bgen_to_plink.%A_%a.out
#SBATCH  --error=run_logs/run_bgen_to_plink.%A_%a.err
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=100000
#SBATCH --time=1-00:00:00
#SBATCH -p mrivas,normal

set -beEuo pipefail

# job start header (for use with array-job module)
_SLURM_JOBID=${SLURM_JOBID:=0} # use 0 for default value (for debugging purpose)
_SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:=1}
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname = $(hostname) SLURM_JOBID = ${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${_SLURM_ARRAY_TASK_ID}" >&2

src="/oak/stanford/groups/mrivas/ukbb24983/imp/pgen/bgen2pgen.sh"

bash ${src} ${_SLURM_ARRAY_TASK_ID} 8 100000

# job finish footer (for use with array-job module)
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname = $(hostname) SLURM_JOBID = ${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${_SLURM_ARRAY_TASK_ID}" >&2


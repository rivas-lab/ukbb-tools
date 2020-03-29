#!/bin/bash
#SBATCH --job-name=bgen2pgen
#SBATCH --output=logs/1_bgen2pgen.%A_%a.out
#SBATCH  --error=logs/1_bgen2pgen.%A_%a.err
#SBATCH --nodes=1
#SBATCH --cores=4
#SBATCH --mem=40000
#SBATCH --time=1:00:00
#SBATCH -p mrivas,normal

set -beEuo pipefail

# job start header (for use with array-job module)
_SLURM_JOBID=${SLURM_JOBID:=0} # use 0 for default value (for debugging purpose)
_SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:=1}
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname = $(hostname) SLURM_JOBID = ${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${_SLURM_ARRAY_TASK_ID}" >&2

src="1_hap_bgen2pgen.sh"

if [ ${_SLURM_ARRAY_TASK_ID} -eq 23 ] ; then
    chr="X"
elif [ ${_SLURM_ARRAY_TASK_ID} -eq 24 ] ; then
    chr="XY"
else
    chr="${_SLURM_ARRAY_TASK_ID}"
fi

bash ${src} ${chr} 4 30000

# job finish footer (for use with array-job module)
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname = $(hostname) SLURM_JOBID = ${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${_SLURM_ARRAY_TASK_ID}" >&2


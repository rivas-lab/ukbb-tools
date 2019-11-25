#!/bin/bash
#SBATCH    --job-name=GRM.combine
#SBATCH --output=logs/GRM.combine.%A.out
#SBATCH  --error=logs/GRM.combine.%A.err
#SBATCH  --nodes=1
#SBATCH  --cores=1
#SBATCH    --mem=8000
#SBATCH   --time=12:00:00
#SBATCH -p mrivas,normal

set -beEuo pipefail

############################################################
# config params
############################################################

cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
source "/oak/stanford/groups/mrivas/users/$USER/repos/rivas-lab/ukbb-tools/15_GRM/15_GRM_misc.sh"

############################################################
# job start
############################################################

_SLURM_JOBID=${SLURM_JOBID:=0} # use 0 for default value (for debugging purpose)
_SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:=1}
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname=$(hostname) ; SLURM_JOBID=${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID=${_SLURM_ARRAY_TASK_ID}" >&2

############################################################
# body
############################################################

pop=white_british
gcta_grm_combine ${pop}

############################################################
# job finish footer (for use with array-job module)
############################################################
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname=$(hostname) ; SLURM_JOBID=${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID=${_SLURM_ARRAY_TASK_ID}" >&2


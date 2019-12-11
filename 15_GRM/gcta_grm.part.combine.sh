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

pop="white_british"
n_parts=100
source "/oak/stanford/groups/mrivas/users/$USER/repos/rivas-lab/ukbb-tools/15_GRM/15_GRM_misc.sh"

############################################################
# job start
############################################################

_SLURM_JOBID=${SLURM_JOBID:=0} # use 0 for default value (for debugging purpose)
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname=$(hostname) ; SLURM_JOBID=${_SLURM_JOBID}" >&2

############################################################
# body
############################################################

gcta_grm_combine ${pop} ${n_parts}

############################################################
# job finish footer (for use with array-job module)
############################################################
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname=$(hostname) ; SLURM_JOBID=${_SLURM_JOBID}" >&2


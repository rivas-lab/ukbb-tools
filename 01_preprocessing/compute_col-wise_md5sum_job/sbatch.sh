#!/bin/bash

#SBATCH --job-name=ukb_md5
#SBATCH   --output=ukb_md5.%A_%a.out
#SBATCH    --error=ukb_md5.%A_%a.err
#SBATCH --time=1-0:00:00
#SBATCH --qos=normal
#SBATCH -p owners,normal,mrivas
#SBATCH --nodes=1
#SBATCH --cores=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#################
# Usage: $ sbatch --array=1-2,4%1 sbatch.sh
#
set -beEu -o pipefail

# automatically get cores and mem settings from the above
cores=$( cat $0 | egrep '^#SBATCH --cores=' | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='   | awk -v FS='=' '{print $NF}' )

echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}" >&2

task_id=$SLURM_ARRAY_TASK_ID
bash task.sh $task_id $mem $cores

echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}" >&2


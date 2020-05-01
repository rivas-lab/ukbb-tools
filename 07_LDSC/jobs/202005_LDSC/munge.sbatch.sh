#!/bin/bash
#SBATCH --job-name=munge
#SBATCH --output=logs/munge.%A_%a.out
#SBATCH  --error=logs/munge.%A_%a.err
#SBATCH --nodes=1
#SBATCH --mem=6000
#SBATCH --cores=1
#SBATCH --time=1-0:00:00

# set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.0.1"
NUM_POS_ARGS="1"

############################################################
# config params
############################################################
src_dir="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/07_LDSC/jobs/202005_LDSC"
icdinfo="${src_dir}/GBE_ID.lst"
batch_size=100
batch_idx=${SLURM_ARRAY_TASK_ID:=1}

get_start_idx () {
    local batch_size=$1
    local batch_idx=$2
    perl -e "print( ( ${batch_idx}  - 1 ) * ${batch_size} + 1 )"
}

get_end_idx () {
    local batch_size=$1
    local batch_idx=$2    
    perl -e "print( ${batch_idx} * ${batch_size} )"
}

get_GBE_IDs () {
    local batch_size=$1
    local batch_idx=$2
    local icdinfo=$3
    
    cat ${icdinfo} | cut -f1 \
    | awk -v s=$(get_start_idx ${batch_size} ${batch_idx}) \
    -v e=$(get_end_idx ${batch_size} ${batch_idx}) \
    's <= NR && NR <= e'
}

############################################################
# job start
############################################################

echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname=$(hostname); SLURM_JOBID=${SLURM_JOBID:=0}; SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:=1}" >&2

############################################################
# body
############################################################

get_GBE_IDs ${batch_size} ${batch_idx} ${icdinfo} | while read GBE_ID ; do
    echo ${GBE_ID}
    bash ${src_dir}/munge.sh ${GBE_ID}
done

############################################################
# job finish footer (for use with array-job module)
############################################################
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname=$(hostname); SLURM_JOBID=${SLURM_JOBID}; SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2

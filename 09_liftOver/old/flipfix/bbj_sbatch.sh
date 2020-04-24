#!/bin/bash

#SBATCH --job-name=flipfix
#SBATCH   --output=logs/flipfix.%A_%a.out
#SBATCH    --error=logs/flipfix.%A_%a.err
#SBATCH --time=0:30:00
#SBATCH --qos=normal
#SBATCH -p owners,normal,mrivas
#SBATCH --nodes=1
#SBATCH --cores=1
#SBATCH --mem=6000
#SBATCH --mail-type=END,FAIL
#################
set -beEu -o pipefail
_SLURM_JOBID=${SLURM_JOBID:=0} # use 0 for default value (for debugging purpose)
_SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:=1}
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname = $(hostname) SLURM_JOBID = ${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${_SLURM_ARRAY_TASK_ID}" >&2
task_id=${_SLURM_ARRAY_TASK_ID}
#################
# set constants

input_file_list="/home/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/flipfix/jobs/bbj/input.files.lst"

ref_fa="/oak/stanford/groups/mrivas/public_data/genomes/hg19/hg19.fa"

src="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/flipfix/flipfix.sh"

#################
# set filenames

in_file=$(cat $input_file_list | awk -v nr=$task_id 'NR==nr')
out_file=$(echo ${in_file} | sed -e "s%bbj/combined_aut_X%bbj/flipfixed_combined_aut_X%g")

if [ ! -d $(dirname $out_file) ] ; then mkdir -p $(dirname $out_file) ; fi
bash $src $in_file $ref_fa | bgzip -l 9 > ${out_file%.gz}.gz

#################
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname = $(hostname) SLURM_JOBID = ${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${_SLURM_ARRAY_TASK_ID}" >&2



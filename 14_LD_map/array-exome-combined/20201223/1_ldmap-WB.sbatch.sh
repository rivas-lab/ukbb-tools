#!/bin/bash
#SBATCH --job-name=LD_MAP
#SBATCH --output=logs/ldmap.%A_%a.out
#SBATCH  --error=logs/ldmap.%A_%a.err
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=40000
#SBATCH --time=7-00:00:00
#SBATCH -p mrivas
#SBATCH --qos high_p

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.1.0"

set -beEuo pipefail
############################################################

source /oak/stanford/groups/mrivas/users/${USER}/repos/rivas-lab/ukbb-tools/14_LD_map/14_LD_map_misc.sh
ml load plink plink2

cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )

# job start header (for use with array-job module)
_SLURM_JOBID=${SLURM_JOBID:=0} # use 0 for default value (for debugging purpose)
_SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:=24}
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname=$(hostname) ; SLURM_JOBID=${_SLURM_JOBID}" >&2
chr_idx=${_SLURM_ARRAY_TASK_ID}

pop="white_british"

#########################
chr=$(echo ${chr_idx} | sed -e 's/23/X/g' | sed -e 's/24/Y/g' | sed -e 's/25/XY/g' | sed -e 's/26/MT/g')
#########################
memory=$(perl -e "print(int(${mem} * 0.8))")
# for plink, we allocate the 80% of the memory
bfile="/scratch/groups/mrivas/ukbb24983/array-exome-combined/pgen/ukb24983_cal_hla_cnv_exomeOQFE"
ldmap_d="/scratch/groups/mrivas/ukbb24983/array-exome-combined/ldmap/ldmap_20201223"
out_prefix="${ldmap_d}/ukb24983_cal_hla_exomeOQFE.${pop}.chr${chr}"
keep_file="/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200828/ukb24983_${pop}.phe"
# plink_opts="--extract:${ldmap_d}/ukb24983_cal_hla_imp.variants.lst" # one can put mat threshold etc.
plink_opts="--chr:${chr}"

# compute_ld_map_wrapper \
compute_ld_r2_wrapper \
    "${bfile}" \
    "${out_prefix}" \
    "${cores}" "${memory}" \
    "${keep_file}" \
    "${plink_opts}"

# job finish footer (for use with array-job module)
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname=$(hostname) ; SLURM_JOBID=${_SLURM_JOBID}" >&2

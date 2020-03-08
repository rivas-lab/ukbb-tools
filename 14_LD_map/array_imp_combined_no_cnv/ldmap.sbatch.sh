#!/bin/bash
#SBATCH --job-name=LD_MAP
#SBATCH --output=logs/ldmap.%A_%a.out
#SBATCH  --error=logs/ldmap.%A_%a.err
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=64000
#SBATCH --time=2-00:00:00
#SBATCH -p mrivas,normal

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.0.1"

set -beEuo pipefail
############################################################

source /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/14_LD_map/14_LD_map_misc.sh
#source "$(dirname ${SRCDIR})/14_LD_map_misc.sh"
ml load plink plink2

cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )

# job start header (for use with array-job module)
_SLURM_JOBID=${SLURM_JOBID:=0} # use 0 for default value (for debugging purpose)
_SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:=1}
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname=$(hostname) ; SLURM_JOBID=${_SLURM_JOBID}" >&2
pop_idx=${_SLURM_ARRAY_TASK_ID}

pops=(
"white_british"
"non_british_white"
"african"
"s_asian"
"e_asian"
)

#########################
pop=$(echo ${pops[@]} | awk -v idx=${pop_idx} '{print $idx}')
#pop="white_british"
#pop="e_asian"
#########################
memory=$(perl -e "print(int(${mem} * 0.8))")
# for plink, we allocate the 80% of the memory
bfile="/scratch/groups/mrivas/ukbb24983/array_imp_combined_no_cnv/pgen/ukb24983_ukb24983_cal_hla_imp"
# bfile="/oak/stanford/groups/mrivas/ukbb24983/array_imp_combined_no_cnv/pgen/ukb24983_ukb24983_cal_hla_imp"
out_prefix="/scratch/groups/mrivas/ukbb24983/array_imp_combined_no_cnv/ldmap/ukb24983_cal_hla_imp.${pop}"
keep_file="/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_${pop}.phe"
plink_opts="" # one can put mat threshold etc.

compute_ld_map_wrapper \
    "${bfile}" \
    "${out_prefix}" \
    "${cores}" "${memory}" \
    "${keep_file}" \
    "${plink_opts}"

# job finish footer (for use with array-job module)
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname=$(hostname) ; SLURM_JOBID=${_SLURM_JOBID}" >&2

#!/bin/bash
#SBATCH --job-name=LD_MAP
#SBATCH --output=logs/ldmap.%A_%a.out
#SBATCH  --error=logs/ldmap.%A_%a.err
#SBATCH --nodes=1
#SBATCH --cores=6
#SBATCH --mem=48000
#SBATCH --time=2-00:00:00
#SBATCH -p mrivas,normal

set -beEuo pipefail

cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )

plink_opts="--memory $(perl -e "print(int(${mem} * 0.8))") --threads ${cores}"

#########################
ml load plink plink2
pop=$1
#pop="white_british"
#pop="e_asian"
#########################

compute_ld_map () {
    local chr=$1
    local pop=$2
    
    local keep_file="${OAK}/ukbb24983/sqc/population_stratification/ukb24983_${pop}.phe"
    
    local dataset_dir="${OAK}/ukbb24983/imp"
    local bfile="${dataset_dir}/pgen/ukb24983_imp_chr${chr}_v3"
    local out_prefix="${dataset_dir}/ldmap_common/ukb24983_imp_common_chr${chr}_v3.${pop}"

    if [ ! -d $(dirname ${out_prefix}) ] ; then mkdir -p $(dirname ${out_prefix}) ; fi

    plink2 ${plink_opts} \
        --maf 0.01 \
        --pfile ${bfile} vzs --keep ${keep_file} \
        --allow-extra-chr \
        --indep-pairwise 50 5 0.5 \
        --out ${out_prefix}.bool

    plink ${plink_opts} \
        --maf 0.01 \
        --bfile ${bfile} --keep ${keep_file} \
        --allow-extra-chr \
        --ld-window-kb 1000 --ld-window-r2 0.1 --r2 gz \
        --out ${out_prefix}.ld_map    
}

# job start header (for use with array-job module)
_SLURM_JOBID=${SLURM_JOBID:=0} # use 0 for default value (for debugging purpose)
_SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:=22}
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname=$(hostname) ; SLURM_JOBID=${_SLURM_JOBID}" >&2
chr=${_SLURM_ARRAY_TASK_ID}

#########################
compute_ld_map ${chr} ${pop}
#########################

# job finish footer (for use with array-job module)
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname=$(hostname) ; SLURM_JOBID=${_SLURM_JOBID}" >&2

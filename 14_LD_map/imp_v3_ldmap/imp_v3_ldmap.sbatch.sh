#!/bin/bash
#SBATCH --job-name=LD_MAP
#SBATCH --output=imp_v3_ldmap_logs/ldmap.%A_%a.out
#SBATCH  --error=imp_v3_ldmap_logs/ldmap.%A_%a.err
#SBATCH --nodes=1
#SBATCH --cores=6
#SBATCH --mem=51200
#SBATCH --time=2-00:00:00
#SBATCH -p mrivas,normal

set -beEuo pipefail

#########################
pop=$1
#pop="white_british"
#pop="e_asian"
#########################

array_job_idx_to_chrom () {
    local idx=$1

    if [ $idx -lt 1 ] || [ $idx -gt 24 ] ; then
        echo "The specified index ($idx) is not supported (it should be 1-24, where 23 = X and 24 = XY)" >&2
        exit 1
    elif [ $idx -eq 23 ] ; then 
        echo "X" 
    elif [ $idx -eq 24 ] ; then
        echo "XY"
    else
        echo "$idx"
    fi
}

compute_ld_map () {
    local chrom=$1
    local pop=$2
    
    local out_dir="/scratch/groups/mrivas/ukbb/24983/imp/ldmap"
    local ukbb_dir="$OAK/ukbb24983"
    local keep_file="${ukbb_dir}/sqc/population_stratification/ukb24983_${pop}.phe"
    local bfile="${ukbb_dir}/imp/pgen/var_QC/ukb24983_imp_chr${chrom}_v3_QC"
    local out_prefix="${out_dir}/ukb24983_imp_chr${chrom}_v3_QC.${pop}"

    ml load plink

    if [ ! -d ${out_dir} ] ; then mkdir -p ${out_dir} ; fi

    plink --allow-extra-chr --bfile ${bfile} --keep ${keep_file} \
        --ld-window-kb 1000 --ld-window-r2 0.1 --r2 gz \
        --out ${out_prefix}.ld_map
    
    plink --allow-extra-chr --bfile ${bfile} --keep ${keep_file} \
        --indep 50 5 2 \
        --out ${out_prefix}.bool 
}

# job start header (for use with array-job module)
_SLURM_JOBID=${SLURM_JOBID:=0} # use 0 for default value (for debugging purpose)
_SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:=24}
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname = $(hostname) ; pop = ${pop} ; SLURM_JOBID = ${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${_SLURM_ARRAY_TASK_ID}" >&2

#########################
chrom=$(array_job_idx_to_chrom ${_SLURM_ARRAY_TASK_ID})
compute_ld_map ${chrom} ${pop}
#########################

# job finish footer (for use with array-job module)
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname = $(hostname) ; pop = ${pop} ; SLURM_JOBID = ${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${_SLURM_ARRAY_TASK_ID}" >&2

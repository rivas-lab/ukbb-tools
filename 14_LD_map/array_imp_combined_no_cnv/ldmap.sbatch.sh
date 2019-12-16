#!/bin/bash
#SBATCH --job-name=LD_MAP
#SBATCH --output=logs/ldmap.%A.out
#SBATCH  --error=logs/ldmap.%A.err
#SBATCH --nodes=1
#SBATCH --cores=12
#SBATCH --mem=192000
#SBATCH --time=2-00:00:00
#SBATCH -p mrivas,normal

set -beEuo pipefail

cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )

out="/oak/stanford/groups/mrivas/ukbb/24983/array_imp_combined_no_cnv/pgen/ukb24983_ukb24983_cal_hla_imp"

plink_opts="--memory $(perl -e "print(int(${mem} * 0.8))") --threads ${cores}"


#########################
ml load plink plink2
pop=$1
#pop="white_british"
#pop="e_asian"
#########################

compute_ld_map () {
    local pop=$1

    local keep_file="${OAK}/ukbb24983/sqc/population_stratification/ukb24983_${pop}.phe"
    
    local dataset_dir="${OAK}/ukbb24983/array_imp_combined_no_cnv"
    local bfile="${dataset_dir}/pgen/ukb24983_ukb24983_cal_hla_imp"
    local out_prefix="${dataset_dir}/ldmap/ukb24983_ukb24983_cal_hla_imp.${pop}"

    if [ ! -d $(dirname ${out_prefix}) ] ; then mkdir -p $(dirname ${out_prefix}) ; fi

    plink2 ${plink_opts} \
        --pfile ${bfile} vzs --keep ${keep_file} \
        --allow-extra-chr \
        --indep-pairwise 50 5 0.5 \
        --out ${out_prefix}.bool

    plink ${plink_opts} \
        --bfile ${bfile} --keep ${keep_file} \
        --allow-extra-chr \
        --ld-window-kb 1000 --ld-window-r2 0.1 --r2 gz \
        --out ${out_prefix}.ld_map    
}

# job start header (for use with array-job module)
_SLURM_JOBID=${SLURM_JOBID:=0} # use 0 for default value (for debugging purpose)
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname=$(hostname) ; SLURM_JOBID=${_SLURM_JOBID}" >&2

#########################
compute_ld_map ${pop}
#########################

# job finish footer (for use with array-job module)
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname=$(hostname) ; SLURM_JOBID=${_SLURM_JOBID}" >&2

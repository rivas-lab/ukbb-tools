#!/bin/bash
#SBATCH --job-name=GWAS_EXOME
#SBATCH --output=rerun_logs/run_exome.%A_%a.out
#SBATCH  --error=rerun_logs/run_exome.%A_%a.err
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=51200
#SBATCH --time=2-00:00:00
#SBATCH -p normal,owners,mrivas

set -beEuo pipefail

# define functions
usage () {
    echo "$0: GWAS re-run script for the exome genotype data"
    echo "usage: sbatch --array=1-<number of array jobs> $0 start_idx (inclusive.)"
}

software_versions () {
    which plink2
    which bgzip
    which python
}

get_field_from_pop () {
    local pop=$1
    local const=7
    echo "white_british non_british_white african e_asian s_asian related others" \
        | tr " " "\n" | awk -v pop=$pop -v const=$const '($0 == pop){print NR + const}'
}

find_phe_path () {
    local info_file=$1
    local min_N=$2
    local col=$3
    local start_idx=$4
    local this_idx=$5

    cat $info_file | awk -v min_N="${min_N}" -v col=$col 'NR > 1 && $col >= min_N' \
        | egrep -v MED \
        | awk -v start_idx=$start_idx -v this_idx=$this_idx 'NR==(start_idx + this_idx - 1) {print $NF}'
}

# get core and memory settings from the header -- passed to gwas script below
cores=$(              cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(                cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )
log_dir=$( dirname $( cat $0 | egrep '^#SBATCH --output=' | awk -v FS='=' '{print $NF}' ))

# check number of command line args and dump usage if that's not right
if [ $# -lt 1 ] ; then usage >&2 ; exit 1 ; fi

# start index
start_idx=$1

# population
if [ $# -gt 1 ] ; then pop=$2 ; else pop="white_british" ; fi
if [ ${pop} != "white_british" ] && [ ${pop} != "non_british_white" ] && [ ${pop} != "african" ] && [ ${pop} != "e_asian" ] && [ ${pop} != "s_asian" ] && [ ${pop} != "related" ] && [ ${pop} != "others" ] ; then
    echo "unsupported population: ${pop}" >&2 ; exit 1 
fi
field=$(get_field_from_pop $pop)

# load sofware, dump which versions are used
export MODULEPATH="/home/groups/mrivas/.modules:${MODULEPATH}"
ml load htslib
ml load python/3.6.1

if grep -q "CPU_GEN:HSW\|CPU_GEN:BDW\|CPU_GEN:SKX" <(a=$(hostname); sinfo -N -n ${a::-4} --format "%50f"); then
   # AVX2 is suitable for use on this node if CPU is recent enough
   ml load plink2/20200727
else
   ml load plink2/20200727-non-AVX2
fi

software_versions >&2

# job start header (for use with array-job module)
_SLURM_JOBID=${SLURM_JOBID:=0} # use 0 for default value (for debugging purpose)
_SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:=1}
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname = $(hostname) SLURM_JOBID = ${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${_SLURM_ARRAY_TASK_ID} ; pop=${pop}" >&2

# get phenotypes to run

min_N_count=10
phenotype_info_file="../05_gbe/exome_phenotype_info.tsv"

echo $phenotype_info_file $min_N_count $field $start_idx $_SLURM_ARRAY_TASK_ID
phe_path=$(find_phe_path ${phenotype_info_file} ${min_N_count} ${field} ${start_idx} ${_SLURM_ARRAY_TASK_ID})
gbeId=$(basename ${phe_path} .phe)

# run exome gwas with default GBE parameters
gwas_out_dir=$(echo $(dirname $(dirname $phe_path)) | awk '{gsub("phenotypedata","exome/gwas"); print}')/${pop}
symlink_dir="/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/current/${pop}"

if [ ! -d ${gwas_out_dir}/logs ] ; then mkdir -p ${gwas_out_dir}/logs ; fi
if [ ! -d $log_dir ] ; then mkdir -p $log_dir ; fi

/share/software/user/open/python/3.6.1/bin/python3 gwas.py --run-exome --run-now --memory $mem --cores $cores --pheno $phe_path --out $gwas_out_dir --population $pop --log-dir $log_dir

# introduce symlinks
file_prefix=ukb24983_v2_hg38.${gbeId}.exome-spb
for ending in "logistic.hybrid" "linear"; do
    if [ -f ${gwas_out_dir}/${file_prefix}.glm.${ending}.gz ]; then
        ln -sf ${gwas_out_dir}/${file_prefix}.glm.${ending}.gz ${symlink_dir}/${file_prefix}.glm.${ending}.gz
    fi
done

ln -sf ${gwas_out_dir}/logs/${file_prefix}.log ${symlink_dir}/logs/${file_prefix}.log

# job finish footer (for use with array-job module)
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname = $(hostname) SLURM_JOBID = ${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${_SLURM_ARRAY_TASK_ID} ; pop=${pop}" >&2

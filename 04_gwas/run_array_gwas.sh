#!/bin/bash
#SBATCH --job-name=RL_ARRAY
#SBATCH --output=rerun_logs/run_array.%A_%a.out
#SBATCH  --error=rerun_logs/run_array.%A_%a.err
#SBATCH --mem=64000
#SBATCH --cores=4
#SBATCH --time=2-00:00:00
#SBATCH -p normal,owners
#SBATCH --constraint=CPU_GEN:HSW|CPU_GEN:BDW|CPU_GEN:SKX, # plink2 avx2 compatibility

set -beEuo pipefail

# define functions
usage () {
    echo "$0: GWAS re-run script for the array genotype data"
    echo "usage: sbatch --array=1-<number of array jobs> $0 start_idx"
    echo ''
    echo "  Before submitting the job, please load htslib and plink2 modules:"
    echo '  $ export MODULEPATH="$HOME/.modules:/home/groups/mrivas/.modules:$OAK/.modules:$MODULEPATH"'
    echo '  $ ml load htslib'
    echo '  $ ml load plink2'    
    echo ''
    echo '  You may check the status of the job (which jobs are finished) using the array-job module:'
    echo '  $ ml load array-job'
    echo '  $ array-job-find_ok.sh rerun_logs'
}

software_versions () {
    which plink2
    which bgzip
}

# automatically get cores, mem, and time settings from the above
batch_cores=$( cat $0 | egrep '^#SBATCH --cores=' | awk -v FS='=' '{print $NF}' )
batch_mem=$(   cat $0 | egrep '^#SBATCH --mem='   | awk -v FS='=' '{print $NF}' )
batch_time=$(  cat $0 | egrep '^#SBATCH --time='  | awk -v FS='=' '{print $NF}' )

# check number of command line args and dump usage if that's not right
if [ $# -lt 1 ] ; then usage >&2 ; exit 1 ; fi

# dump which softwares are used
software_versions >&2

echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}" >&2

# get phenotypes to run
start_idx=$1
this_idx=$SLURM_ARRAY_TASK_ID

phe_path=$(awk -v a=$start_idx '(a <= NR){print $NF}' ../05_gbe/phenotype_info.tsv | awk -v nr=$this_idx 'NR==nr')
gbeId=$(basename $phe_path | awk '{gsub(".phe","");print}')

# run array gwas with default GBE parameters
pop="white_british"
gwasOutDir=$(echo $(dirname $(dirname $phe_path)) | awk '{gsub("phenotypedata","cal/gwas"); print}')/${pop}
if [ ! -d ${gwasOutDir}/logs ] ; then mkdir -p ${gwasOutDir}/logs ; fi
if [ ! -d $(dirname $(readlink -f $0))/rerun_logs ] ; then mkdir -p $(dirname $(readlink -f $0))/rerun_logs ; fi

python gwas.py --run-array --run-now --batch-memory ${batch_mem} --batch-time ${batch_time} --batch-cores ${batch_cores} --pheno $phe_path --out $gwasOutDir --population $pop --log-dir rerun_logs

# move log file and bgzip output
for type in genotyped; do 
    for ending in "logistic.hybrid" "linear"; do
        if [ -f ${gwasOutDir}/ukb24983_v2_1.${gbeId}.${type}.glm.${ending} ]; then
            bgzip -f ${gwasOutDir}/ukb24983_v2_1.${gbeId}.${type}.glm.${ending}
        fi
    done
    if [ -f ${gwasOutDir}/ukb24983_v2_1.${gbeId}.${type}.log ]; then
        mv ${gwasOutDir}/ukb24983_v2_1.${gbeId}.${type}.log ${gwasOutDir}/logs/
    fi    
done

echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}" >&2

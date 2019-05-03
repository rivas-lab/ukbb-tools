#!/bin/bash
#SBATCH --job-name=RL_ARRAY
#SBATCH --output=rerun_logs/run_array.%A_%a.out
#SBATCH  --error=rerun_logs/run_array.%A_%a.err
#SBATCH --mem=24000
#SBATCH --cores=4
#SBATCH --time=2-00:00:00
#SBATCH -p normal,owners

set -beEuo pipefail

# define functions
usage () {
    echo "$0: GWAS re-run script for the array genotype data"
    echo "usage: sbatch --array=1-<number of array jobs> $0 start_idx"
    echo ''
    echo '  You may check the status of the job (which jobs are finished) using the array-job module:'
    echo '  $ ml load array-job'
    echo '  $ array-job-find_ok.sh rerun_logs'
}

software_versions () {
    which plink2
    which bgzip
}

# get core and memory settings from the header -- passed to gwas script below
cores=$(              cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(                cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )
log_dir=$( dirname $( cat $0 | egrep '^#SBATCH --output=' | awk -v FS='=' '{print $NF}' ))

# check number of command line args and dump usage if that's not right
if [ $# -lt 1 ] ; then usage >&2 ; exit 1 ; fi

# load sofware, dump which versions are used
export MODULEPATH="/home/groups/mrivas/.modules:${MODULEPATH}"
ml load htslib

if grep -q "CPU_GEN:HSW\|CPU_GEN:BDW\|CPU_GEN:SKX" <(a=$(hostname); sinfo -N -n ${a::-4} --format "%50f"); then
   # AVX2 is suitable for use on this node if CPU is recent enough
   ml load plink2/20190402
else
   ml load plink2/20190402-non-AVX2
fi

software_versions >&2

# job start header (for use with array-job module)
_SLURM_JOBID=${SLURM_JOBID:=0}
_SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:=1}
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname = $(hostname) SLURM_JOBID = ${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${_SLURM_ARRAY_TASK_ID}" >&2

# get phenotypes to run
start_idx=$1
this_idx=$_SLURM_ARRAY_TASK_ID

phe_path=$(awk -v a=$start_idx '(a <= NR){print $NF}' ../05_gbe/phenotype_info.tsv | awk -v nr=$this_idx 'NR==nr')
gbeId=$(basename $phe_path | awk '{gsub(".phe","");print}')

# run array gwas with default GBE parameters
pop="white_british"
gwasOutDir=$(echo $(dirname $(dirname $phe_path)) | awk '{gsub("phenotypedata","cal/gwas"); print}')/${pop}
if [ ! -d ${gwasOutDir}/logs ] ; then mkdir -p ${gwasOutDir}/logs ; fi
if [ ! -d $log_dir ] ; then mkdir -p $log_dir ; fi

python gwas.py --run-array --run-now --memory $mem --cores $cores --pheno $phe_path --out $gwasOutDir --population $pop --log-dir $log_dir

# move log file and bgzip output
for type in genotyped; do 
    for ending in "logistic.hybrid" "linear"; do
        if [ -f ${gwasOutDir}/ukb24983_v2.${gbeId}.${type}.glm.${ending} ]; then
            bgzip --compress-level 9 -f ${gwasOutDir}/ukb24983_v2.${gbeId}.${type}.glm.${ending}
        fi
    done
    if [ -f ${gwasOutDir}/ukb24983_v2.${gbeId}.${type}.log ]; then
        mv -f ${gwasOutDir}/ukb24983_v2.${gbeId}.${type}.log ${gwasOutDir}/logs/
    fi    
done

# job finish footer (for use with array-job module)
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname = $(hostname) SLURM_JOBID = ${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${_SLURM_ARRAY_TASK_ID}" >&2


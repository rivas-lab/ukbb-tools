#!/bin/bash
set -beEuo pipefail

mem=16000
cores=4

pop=$1 # related | others
batch_idx=${SLURM_ARRAY_TASK_ID:=1}

job_idx="20200627-gwas-additional-pops-redo.job.index.tsv"
pheno="/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.20200522.phe"

echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname=$(hostname); SLURM_JOBID=${SLURM_JOBID:=0}; SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:=1}" >&2

################################################
# load the modules
################################################

# export MODULEPATH="$HOME/.modules:/labs/mrivas/.modules:$MODULEPATH"
is_scg=$(echo scg $(hostname) | tr " " "\n" | grep scg | wc -l)
avx2_flag=$( cat /proc/cpuinfo  | grep flags | uniq  | awk -v FS=':' '{print $2}' | tr " " "\n" | cat /dev/stdin <(echo avx2) | grep -i avx2 | wc -l)
if [ "${is_scg}" -eq 2 ] ; then
    export MODULEPATH="$HOME/.modules:/labs/mrivas/.modules:$MODULEPATH"
    ml load htslib python/3.8.2
else
    ml load htslib python/3.6.1
fi
if [ "${avx2_flag}" -eq 2 ] ; then
    ml load plink2/20200314
else
    ml load plink2/20200314-non-AVX2
fi

################################################
# functions
################################################

get_plink_suffix () {
    local GBE_ID=$1

    GBE_CAT=$(echo $GBE_ID | sed -e "s/[0-9]//g")

    if [ "${GBE_CAT}" == "QT_FC" ] || [ "${GBE_CAT}" == "INI" ] ; then
        echo glm.linear
    else
        echo glm.logistic.hybrid
    fi
}

################################################
# configure the params
################################################

if   [ "${pop}" == "related" ] ; then
    out="/oak/stanford/groups/mrivas/projects/related/"
    keep="/oak/stanford/groups/mrivas/ukbb24983/sqc/relatedness_20200514/semi_related.fam"
elif [ "${pop}" == "others" ] ; then
    out="/oak/stanford/groups/mrivas/projects/gwas_others/"
    keep="/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200522/ukb24983_others.phe"
else
    echo "unsupported population ${pop}" >&2 ; exit 1
fi

GBE_ID=$(    cat ${job_idx} | awk -v pop=${pop} '$1==pop' | awk -v idx=${batch_idx} '(NR==idx){print $2}')
pheno_name=$(cat ${job_idx} | awk -v pop=${pop} '$1==pop' | awk -v idx=${batch_idx} '(NR==idx){print $3}')
plink_suffix=$(get_plink_suffix ${GBE_ID})

out1=${out}/ukb24983_v2_hg19.${pheno_name}.array-combined.one_array.${GBE_ID}.${plink_suffix}
out2=${out}/ukb24983_v2_hg19.${pheno_name}.array-combined.both_arrays.${GBE_ID}.${plink_suffix}
out_gz=${out}/ukb24983_v2_hg19.${GBE_ID}.array-combined.${plink_suffix}.gz

log1=${out}/ukb24983_v2_hg19.${pheno_name}.array-combined.one_array.log
log2=${out}/ukb24983_v2_hg19.${pheno_name}.array-combined.both_arrays.log
log=${out}/ukb24983_v2_hg19.${GBE_ID}.array-combined.log

################################################
# run GWAS with plink
################################################

python3 /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/gwas.20200606.py \
    --memory ${mem} --cores ${cores} \
    --run-array-combined --run-now \
    --population all --keep-related --include-x --run-now \
    --pheno ${pheno} \
    --pheno-name ${pheno_name} \
    --out ${out} \
    --keep ${keep}

if [ -f ${log1} ] && [ -f ${log2} ] ; then
# combine the log files
    cat <(echo "#$(basename ${log1})") ${log1} <(echo "#$(basename ${log2})") ${log2} > ${log}
    if [ -f ${log} ] ; then rm ${log1} ${log2} ; fi   
    echo ${log}
fi

if [ -f ${out1} ] && [ -f ${out2} ] ; then
# combine the results files
    {
        ! head -n1 ${out1}
        cat ${out1} ${out2} | grep -v '#' | sort --parallel ${cores} -k1,1V -k2,2n -k3,3
    } | bgzip -l9 -@${cores} > ${out_gz}

    if [ -f ${out_gz} ] ; then rm ${out1} ${out2} ; fi
    
    echo ${out_gz}
fi

echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname=$(hostname); SLURM_JOBID=${SLURM_JOBID}; SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2


exit 0
################################################################################################
# the followings are just some notes
################################################################################################
# instructions for Yosuke
- There are 
  - 219 phes for others
  - 112 phes for related

cd /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/extras/202006-GWAS-finish

sbatch --account=mrivas -p nih_s10 --cpus-per-task=4 --ntasks=1 --mem=16000 --time=7-0:00:00 --job-name=gwas-others --output=logs/gwas-others.%A_%a.out --error=logs/gwas-others.%A_%a.err --array=1-219 20200627-gwas-additional-pops-redo.sh others

Submitted batch job 15887953

scontrol update job 15887953 partition=nih_s10,batch

sbatch --account=mrivas -p nih_s10 --cpus-per-task=4 --ntasks=1 --mem=16000 --time=7-0:00:00 --job-name=gwas-related --output=logs/gwas-related.%A_%a.out --error=logs/gwas-related.%A_%a.err --array=1-112 20200627-gwas-additional-pops-redo.sh related

Submitted batch job 15887954

scontrol update job 15887954 partition=nih_s10,batch

## fix bug (run both one array and both arrays)

 1114  sbatch --account=mrivas -p nih_s10 --cpus-per-task=4 --ntasks=1 --mem=16000 --time=7-0:00:00 --job-name=gwas-others --output=logs/gwas-others.%A_%a.out --error=logs/gwas-others.%A_%a.err --array=1-219 20200627-gwas-additional-pops-redo.sh others
 1115  sbatch --account=mrivas -p nih_s10 --cpus-per-task=4 --ntasks=1 --mem=16000 --time=7-0:00:00 --job-name=gwas-related --output=logs/gwas-related.%A_%a.out --error=logs/gwas-related.%A_%a.err --array=1-112 20200627-gwas-additional-pops-redo.sh related
 1116  scontrol update job 15888144 partition=nih_s10,batch
 1117  scontrol update job 15888145 partition=nih_s10,batch
 
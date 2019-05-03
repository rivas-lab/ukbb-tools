#!/bin/bash
#SBATCH --job-name=RL_EXOME
#SBATCH --output=rerun_logs/run_exome.%A-%a.out
#SBATCH --mem=64000
#SBATCH --cores=4
#SBATCH --time=2-00:00:00
#SBATCH -p normal,owners
# #SBATCH --constraint=CPU_GEN:HSW|CPU_GEN:BDW|CPU_GEN:SKX, # plink2 avx2 compatibility

# dependencies
export MODULEPATH="/home/groups/mrivas/.modules/:${MODULEPATH}"
ml load biology; ml load htslib; ml load plink2/20190402-non-AVX2

# get phenotypes to run
start_idx=$1
end_idx=$2
this_idx=$SLURM_ARRAY_TASK_ID

phe_path=$(awk -v a=$start_idx -v b=$end_idx '(a <= NR && NR <= b){print $NF}' ../05_gbe/phenotype_info.tsv | awk -v nr=$this_idx 'NR==nr')
gbeId=$(basename $phe_path | awk '{gsub(".phe","");print}')

# run exome gwas with default GBE parameters
pop="white_british"
gwasOutDir=$(echo $(dirname $(dirname $phe_path)) | awk '{gsub("phenotypedata","exome/gwas"); print}')/${pop}
mkdir -p ${gwasOutDir}/logs
mkdir -p rerun_logs

python gwas.py --run-exome --run-now --pheno $phe_path --out $gwasOutDir --population $pop --log-dir rerun_logs

# move log file and bgzip output
for type in genotyped; do 
    if [ -f ${gwasOutDir}/ukb24983_v2_1.${gbeId}.${type}.log ]; then
        mv ${gwasOutDir}/ukb24983_v2_1.${gbeId}.${type}.log ${gwasOutDir}/logs/
    fi
    for ending in "logistic.hybrid" "linear"; do
        if [ -f ${gwasOutDir}/ukb24983_v2_1.${gbeId}.${type}.glm.${ending} ]; then
            bgzip -f ${gwasOutDir}/ukb24983_v2_1.${gbeId}.${type}.glm.${ending}
        fi
    done
done
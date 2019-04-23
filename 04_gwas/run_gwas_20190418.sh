#!/bin/bash
#SBATCH --job-name=RL_GWAS
#SBATCH --output=rerun_logs/rl-gwas.20190405.magu.%A-%a.out)
#SBATCH --mem=16000
#SBATCH --cores=4
#SBATCH --time=2-00:00:00
#SBATCH -p normal,owners
#SBATCH --constraint=CPU_GEN:HSW|CPU_GEN:BDW|CPU_GEN:SKX, # plink2 avx2 compatibility

# dependencies
ml load htslib; ml load plink2/20190402-non-AVX2

# get phenotypes to run
start_idx=$1
end_idx=$2
this_idx=$SLURM_ARRAY_TASK_ID

phe_path=$(awk -v a=$start_idx -v b=$end_idx '(a <= NR && NR <= b){print $NF}' ../05_gbe/phenotype_info.tsv | awk -v nr=$this_idx 'NR==nr')
gbeId=$(basename $phe_path | awk '{gsub(".phe","");print}')

# run gwas with default GBE parameters
pop="white_british"
gwasOutDir=$(echo $(dirname $(dirname $phe_path)) | awk '{gsub("phenotypedata","cal/gwas"); print}')/${pop}
mkdir -p ${gwasOutDir}/logs
mkdir -p rerun_logs

python gwas.py --run-array --run-now --pheno $phe_path --out $gwasOutDir --population $pop --log-dir rerun_logs

# move log file and bgzip output
for type in genotyped; do 
    if [ -f ${gwasOutDir}/ukb24983_v2.${gbeId}.${type}.log ]; then
        mv ${gwasOutDir}/ukb24983_v2.${gbeId}.${type}.log ${gwasOutDir}/logs/
    fi
    for ending in "logistic.hybrid" "linear"; do
        if [ -f ${gwasOutDir}/ukb24983_v2.${gbeId}.${type}.glm.${ending} ]; then
            bgzip ${gwasOutDir}/ukb24983_v2.${gbeId}.${type}.glm.${ending}
        fi
    done
done

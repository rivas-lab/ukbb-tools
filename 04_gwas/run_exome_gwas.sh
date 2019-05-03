#!/bin/bash
#SBATCH --job-name=RL_EXOME
#SBATCH --output=rerun_logs/run_exome.%A-%a.out
#SBATCH --mem=32000
#SBATCH --cores=4
#SBATCH --time=2-00:00:00
#SBATCH -p normal,owners

# software dependencies
export MODULEPATH="/home/groups/mrivas/.modules:${MODULEPATH}"
ml load htslib

if grep -q "CPU_GEN:HSW\|CPU_GEN:BDW\|CPU_GEN:SKX" <(a=$(hostname); sinfo -N -n ${a::-4} --format "%50f"); then
   # AVX2 is suitable for use on this node if CPU is recent enough
   ml load plink2/20190402
else
   ml load plink2/20190402-non-AVX2
fi

# get phenotypes to run
start_idx=$1
this_idx=$SLURM_ARRAY_TASK_ID

phe_path=$(awk -v a=$start_idx '(a <= NR){print $NF}' ../05_gbe/phenotype_info.tsv | awk -v nr=$this_idx 'NR==nr')
gbeId=$(basename $phe_path | awk '{gsub(".phe","");print}')

# run exome gwas with default GBE parameters
pop="white_british"
gwasOutDir=$(echo $(dirname $(dirname $phe_path)) | awk '{gsub("phenotypedata","exome/gwas"); print}')/${pop}
mkdir -p ${gwasOutDir}/logs
mkdir -p rerun_logs

# and core/memory constraints from above
cores=$( cat $0 | egrep '^#SBATCH --cores=' | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='   | awk -v FS='=' '{print $NF}' )

python gwas.py --run-exome --run-now --pheno $phe_path --out $gwasOutDir --population $pop --log-dir rerun_logs --cores $cores --memory $mem

# move log file and bgzip output
for type in genotyped; do 
    if [ -f ${gwasOutDir}/ukb24983_v2.${gbeId}.${type}.log ]; then
        mv -f ${gwasOutDir}/ukb24983_v2.${gbeId}.${type}.log ${gwasOutDir}/logs/
    fi
    for ending in "logistic.hybrid" "linear"; do
        if [ -f ${gwasOutDir}/ukb24983_v2.${gbeId}.${type}.glm.${ending} ]; then
            bgzip -f ${gwasOutDir}/ukb24983_v2.${gbeId}.${type}.glm.${ending}
        fi
    done
done

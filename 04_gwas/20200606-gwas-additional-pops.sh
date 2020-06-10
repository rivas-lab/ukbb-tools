#!/bin/bash
set -beEuo pipefail

mem=16000
cores=4

pop=$1 # related | others
if [ $# -gt 1 ] ; then
    regression_type=$2 # linear
else
    regression_type="logistic"
fi
batch_idx=${SLURM_ARRAY_TASK_ID:=1}
batch_size_bin=4
batch_size_qt=10

# export MODULEPATH="$HOME/.modules:/labs/mrivas/.modules:$MODULEPATH"
avx2_flag=$( cat /proc/cpuinfo  | grep flags | uniq  | awk -v FS=':' '{print $2}' | tr " " "\n" | cat /dev/stdin <(echo avx2) | grep -i avx2 | wc -l)
if [ "${avx2_flag}" -eq 2 ] ; then
    # ml load htslib python/3.8.2 plink2/20200314
    ml load htslib python/3.6.1 plink2/20200314
else
    # ml load htslib python/3.8.2 plink2/20200314-non-AVX2
    ml load htslib python/3.6.1 plink2/20200314
fi

pheno="/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.20200522.phe"
pheno_covar_col_end=94

if   [ "${pop}" == "related" ] ; then
    out="/oak/stanford/groups/mrivas/projects/related/"
    # out="/oak/stanford/groups/mrivas/projects/dev-related/"
    keep="/oak/stanford/groups/mrivas/ukbb24983/sqc/relatedness_20200514/semi_related.fam"
elif [ "${pop}" == "others" ] ; then
    out="/oak/stanford/groups/mrivas/projects/gwas_others/"
    # out="/oak/stanford/groups/mrivas/projects/dev-gwas_others/"
    keep="/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200522/ukb24983_others.phe"
else
    echo "unsupported population ${pop}" >&2 ; exit 1
fi

if [ "${regression_type}" == "linear" ] ; then
    start_idx=$( perl -e "print(int(((${batch_idx} - 1) * ${batch_size_qt}) + 1))" )
    end_idx=$(   perl -e "print(int(  ${batch_idx}      * ${batch_size_qt}))" )
    pheno_name=$(head -n1 ${pheno} | tr "\t" "\n" | egrep -n 'INI|QT_FC' | awk -v pheno_covar_col_end=${pheno_covar_col_end} -v FS=':' '(NR>pheno_covar_col_end){print $1}' | awk -v start_idx=${start_idx} -v end_idx=${end_idx} 'start_idx <= NR && NR <= end_idx' | tr "\n" "," | rev | cut -c2- | rev)
elif [ "${regression_type}" == "logistic" ] ; then
    start_idx=$( perl -e "print(int(((${batch_idx} - 1) * ${batch_size_bin}) + 1))" )
    end_idx=$(   perl -e "print(int(  ${batch_idx}      * ${batch_size_bin}))" )
    pheno_name=$(head -n1 ${pheno} | tr "\t" "\n" | egrep -n -v 'INI|QT_FC' | awk -v pheno_covar_col_end=${pheno_covar_col_end} -v FS=':' '(NR>pheno_covar_col_end){print $1}' | awk -v start_idx=${start_idx} -v end_idx=${end_idx} 'start_idx <= NR && NR <= end_idx' | tr "\n" "," | rev | cut -c2- | rev)
else
    echo "unsupported regression_type ${regression_type}" >&2 ; exit 1
fi

# python3 /oak/stanford/groups/mrivas/users/guhan/repos/ukbb-tools/04_gwas/gwas.py \
python3 /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/gwas.20200606.py \
    --memory ${mem} --cores ${cores} \
    --run-array-combined --run-now \
    --population all --keep-related --include-x --run-now \
    --pheno ${pheno} \
    --pheno-name ${pheno_name} \
    --out ${out} \
    --keep ${keep}

exit 0
# instructions
- There are 1569 QTs and 1899 binary phenotypes
- 156 array jobs for QTs, 475 array jobs for bins
- Yosuke --> all of QTs and binary 1-200
- Manny --> binary 201-475

## scg4

### Yosuke
cd /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas
for pop in related others ; do 
sbatch --account=mrivas -p nih_s10 --cpus-per-task=4 --ntasks=1 --mem=16000 --time=7-0:00:00 --job-name=gwas-${pop} --output=logs/gwas-${pop}.%A_%a.out --error=logs/gwas-${pop}.%A_%a.err --array=1-156 20200606-gwas-additional-pops.sh ${pop} linear

sbatch --account=mrivas -p nih_s10 --cpus-per-task=4 --ntasks=1 --mem=16000 --time=7-0:00:00 --job-name=gwas-${pop} --output=logs/gwas-${pop}.%A_%a.out --error=logs/gwas-${pop}.%A_%a.err --array=1-200 20200606-gwas-additional-pops.sh ${pop}
done

### Manny
cd /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas
for pop in related others ; do 
sbatch --account=mrivas -p nih_s10 --cpus-per-task=4 --ntasks=1 --mem=16000 --time=7-0:00:00 --job-name=gwas-${pop} --output=logs/gwas-${pop}.%A_%a.out --error=logs/gwas-${pop}.%A_%a.err --array=201-475 20200606-gwas-additional-pops.sh ${pop}
done

## Sherlock

### Yosuke
cd /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas
for pop in related others ; do 
sbatch -p mrivas --qos=high_p --cpus-per-task=4 --ntasks=1 --mem=16000 --time=7-0:00:00 --job-name=gwas-${pop} --output=logs/gwas-${pop}.%A_%a.out --error=logs/gwas-${pop}.%A_%a.err --array=1-156 20200606-gwas-additional-pops.sh ${pop} linear

sbatch -p mrivas --qos=high_p --cpus-per-task=4 --ntasks=1 --mem=16000 --time=7-0:00:00 --job-name=gwas-${pop} --output=logs/gwas-${pop}.%A_%a.out --error=logs/gwas-${pop}.%A_%a.err --array=1-200 20200606-gwas-additional-pops.sh ${pop}
done

# Submitted batch job 1887505
# Submitted batch job 1887506
# Submitted batch job 1887507
# Submitted batch job 1887508

### Manny
cd /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas
for pop in related others ; do 
sbatch -p mrivas --qos=high_p --cpus-per-task=4 --ntasks=1 --mem=16000 --time=7-0:00:00 --job-name=gwas-${pop} --output=logs/gwas-${pop}.%A_%a.out --error=logs/gwas-${pop}.%A_%a.err --array=201-475 20200606-gwas-additional-pops.sh ${pop}
done

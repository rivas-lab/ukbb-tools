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
batch_idx=${SLURM_ARRAY_TASK_ID:=201}
batch_size_bin=4
batch_size_qt=10

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
################################################################################################
# the followings are just some notes
################################################################################################
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

scontrol hold 1887508_[115-200] # gwas-others
scontrol hold 1887506_[97-200]  # gwas-related

pop=related
idx_range=97-200
sbatch --account=mrivas -p nih_s10 --cpus-per-task=4 --ntasks=1 --mem=16000 --time=7-0:00:00 --job-name=gwas-${pop} --output=logs/gwas-${pop}.%A_%a.out --error=logs/gwas-${pop}.%A_%a.err --array=${idx_range} 20200606-gwas-additional-pops.sh ${pop}

pop=others
idx_range=118-200
sbatch --account=mrivas -p nih_s10 --cpus-per-task=4 --ntasks=1 --mem=16000 --time=7-0:00:00 --job-name=gwas-${pop} --output=logs/gwas-${pop}.%A_%a.out --error=logs/gwas-${pop}.%A_%a.err --array=${idx_range} 20200606-gwas-additional-pops.sh ${pop}

### Manny
cd /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas
for pop in related others ; do 
sbatch -p mrivas --qos=high_p --cpus-per-task=4 --ntasks=1 --mem=16000 --time=7-0:00:00 --job-name=gwas-${pop} --output=logs/gwas-${pop}.%A_%a.out --error=logs/gwas-${pop}.%A_%a.err --array=201-475 20200606-gwas-additional-pops.sh ${pop}
done

### the original idea was... write something like thie: 20200606-others.yt.sh
python /oak/stanford/groups/mrivas/users/guhan/repos/ukbb-tools/04_gwas/gwas.py \
    --run-array-combined --run-now \
    --memory $mem --cores $cores \
    --pheno /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.20200522.phe --out /oak/stanford/groups/mrivas/projects/gwas_others/ --population all --keep-related --pheno-name 1829-3562  --include-x --run-now --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200522/ukb24983_others.phe
exit 0
sbatch -p mrivas --qos=high_p --nodes=1 --mem=150000 --cores=16 --time=7-0:00:00 --job-name=others --output=logs/others.%A.out --error=logs/others.%A.err 20200606-others.yt.sh
### but this never worked (too slow)


### 2020/6/16 -- it turned out that the Manny's jobs [201-475] failed due to the mv/post_processing
I modified gwas.20200606.py script so that we only run plink for the both_array variants

for pop in related others ; do 
sbatch --account=mrivas -p nih_s10 --cpus-per-task=4 --ntasks=1 --mem=16000 --time=7-0:00:00 --job-name=gwas-${pop} --output=logs/gwas-${pop}.%A_%a.out --error=logs/gwas-${pop}.%A_%a.err --array=201-475 20200606-gwas-additional-pops.sh ${pop}
done

Submitted batch job 15775587
Submitted batch job 15775588

-bash-4.2$ scontrol update job 15775587 partition=nih_s10,batch
-bash-4.2$ scontrol update job 15775588 partition=nih_s10,batch

# in terms of the GWAS finalization... (2020/6/16 0:11 am)

- For others and related,
  - We are now applying "rename" scripts to Yosuke's portion of files. It seems like they have 1080969 lines which is expected number based on the pvar file
  - We submitted the GWAS jobs for the Manny's portion (idx 201-475), focusing on the both_array variants

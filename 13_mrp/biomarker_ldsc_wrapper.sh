#!/bin/bash
set -beEuo pipefail

IFS=$'\r\n' GLOBIGNORE='*' command eval  'gbe_ids=($(cut -f4 ../sumstat_paths.tsv | tail -n +2))'
pop="metal";
geno_dataset="array-combined";
ldsc_dir="/oak/stanford//groups/mrivas/ukbb24983/${geno_dataset}/ldsc"

for ((i = 0; i < ${#gbe_ids[@]}; i++)); do 
    for ((j = i + 1; j < ${#gbe_ids[@]}; j++)); do 
        GBE_ID_1="${gbe_ids[i]}"
        GBE_ID_2="${gbe_ids[j]}"
        UKB_f1="${ldsc_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID_1}.${geno_dataset}.sumstats.gz"
        UKB_f2="${ldsc_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID_2}.${geno_dataset}.sumstats.gz"
        out_f="${pop}.${GBE_ID_1}.${GBE_ID_2}"
        src="/oak/stanford/groups/mrivas/users/${USER}/repos/ukbb-tools/07_LDSC/helpers/ldsc_rg.sh"
        if [ ! -f ${out_f}.log ] ; then
            sbatch ${src} --scratch ${UKB_f1} ${UKB_f2} ${out_f}
        fi
    done
done

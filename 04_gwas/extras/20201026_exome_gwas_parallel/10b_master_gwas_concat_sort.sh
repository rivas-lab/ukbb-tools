#!/bin/bash
set -beEuo pipefail

pop=$1

ml load R/3.6 gcc

data_d="/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE"
out_f=${data_d}/ukb24983_exomeOQFE.${pop}.p1e-3.tsv.gz

if [ "${pop}" == "metal" ] ; then

    Rscript 10_master_gwas_concat_sort.R \
    ${out_f%.gz} \
    $(dirname ${data_d})/$(basename ${data_d})_filtered_p1e-3/${pop}.tsv.gz \

else

    Rscript 10_master_gwas_concat_sort.R \
    ${out_f%.gz} \
    $(dirname ${data_d})/$(basename ${data_d})_filtered_p1e-3/${pop}.linear.tsv.gz \
    $(dirname ${data_d})/$(basename ${data_d})_filtered_p1e-3/${pop}.logistic.hybrid.tsv.gz

fi

bgzip -l9 -@6 ${out_f%.gz}
tabix -s1 -b2 -e2 -c'#' ${out_f}
echo ${out_f}

exit 0
#######################
bash 10b_master_gwas_concat_sort.sh african

for pop in 'white_british' 'non_british_white' 's_asian' 'e_asian' 'related' 'others' 'metal' ; do
sbatch -p mrivas --qos=high_p --nodes=1 --mem=64000 --cores=6 --time=7:00:00 \
    --job-name=concat_sort.${pop} --output=logs/concat_sort.%A.out --error=logs/concat_sort.%A.err \
    10b_master_gwas_concat_sort.sh ${pop}
done


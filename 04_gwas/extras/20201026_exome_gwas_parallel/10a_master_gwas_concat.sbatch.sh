#!/bin/bash
set -beEuo pipefail

batch_idx=$1

data_d="/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_filtered_p1e-3"
pops=('white_british' 'non_british_white' 'african' 's_asian' 'e_asian' 'related' 'others')

src="10_master_gwas_concat.sh"
nCores=6

if [ ${batch_idx} -gt 15 ] ; then
    echo "specified (${batch_idx}) is too large" >&2
    exit 1

elif [ ${batch_idx} -eq 15 ] ; then
    # metal
    pop=metal
    header_phe="INI50"
    out_f=${data_d}/${pop}.tsv.gz
else
    pop_idx=$( perl -e "print(int((${batch_idx} - 1) / 2))")
    pop=${pops[${pop_idx}]}
    if [ $( perl -e "print(${batch_idx} % 2)") -eq 1 ] ; then
        header_phe="INI50"
        out_f=${data_d}/${pop}.linear.tsv.gz
    else
        header_phe="HC269"
        out_f=${data_d}/${pop}.logistic.hybrid.tsv.gz
    fi
fi

echo ${batch_idx} ${pop} ${header_phe} $(basename ${out_f})
    
bash ${src} ${pop} ${header_phe} | bgzip -@${nCores} > ${out_f}

exit 0
#######################

sbatch -p mrivas --qos=high_p --nodes=1 --mem=24000 --cores=6 --time=7-0:00:00 \
    --job-name=concat --output=logs/concat.%A.out --error=logs/concat.%A.err \
    --array=1-15 ${parallel_sbatch_sh} \
    10a_master_gwas_concat.sbatch.sh ${parallel_idx} 1


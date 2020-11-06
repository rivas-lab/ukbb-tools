#!/bin/bash
set -beEuo pipefail

batch_idx=$1
pop=$2
if [ $# -gt 2 ] ; then cores=$3 ; else cores=1 ; fi
# pop="non_british_white"

oak_data_d="/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_scg"
data_d="/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE"
combine_src="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/extras/20201026_exome_gwas_parallel/4_combine_output.sh"

get_trait_name () {
    idx=$1
#     master_phe=/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.20201002.sorted.phe.zst
    master_phe=/scratch/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.20201002.sorted.phe.zst
    ! zstdcat ${master_phe} \
    | head -n1 \
    | cut -f1994-3562 | tr '\t' '\n' \
    | awk -v idx=$idx 'NR==idx'
}

if [ "${batch_idx}" -eq 1570 ] ; then
    # log file
    template_f=${oak_data_d}/${pop}-batch/logs/ukb24983_exomeOQFE.batch__BATCH__.QT_ALL.glm.log
    combined_f=${data_d}/${pop}/logs/ukb24983_exomeOQFE.QT_ALL.glm.log.gz
    if [ ! -s ${combined_f} ] ; then
        bash ${combine_src} \
        --log \
        --cores ${cores} \
        --n_batch $([ "${pop}" == "white_british" ] && echo "1000" || echo "100") \
        ${template_f} ${combined_f}
    fi
elif [ "${batch_idx}" -lt 1570 ] ; then
    # plink files
    trait=$(get_trait_name ${batch_idx})
    template_f=${oak_data_d}/${pop}-batch/ukb24983_exomeOQFE.batch__BATCH__.${trait}.glm.linear
    combined_f=${data_d}/${pop}/ukb24983_exomeOQFE.${trait}.glm.linear.gz
    if [ ! -s ${combined_f} ] ; then
        bash ${combine_src} \
        --cores ${cores} \
        --n_batch $([ "${pop}" == "white_british" ] && echo "1000" || echo "100") \
        ${template_f} ${combined_f}
    else
        echo ${combined_f}
    fi                                
fi

exit 0
########################################################################
bash 4_combine_output_QT_ALL_from_oak.sh 1 african

#pop=african ; \

for pop in 'african' 's_asian' 'e_asian' 'related' 'others' ; do
sbatch -p mrivas,normal \
    --time=4:0:00 --mem=6000 --nodes=1 --cores=1 \
    --job-name=combine_output_QT_ALL \
    --output=logs_scratch/combine_output_QT_ALL.%A_%a.out \
    --error=logs_scratch/combine_output_QT_ALL.%A_%a.err \
    --array=1-157 ${parallel_sbatch_no_err_check_sh} \
    4_combine_output_QT_ALL_from_oak.sh ${parallel_idx} 10 ${pop}
done


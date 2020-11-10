#!/bin/bash
set -beEuo pipefail

pop=$1
trait=$2

if [ $# -gt 2 ] ; then cores=$3 ; else cores=1 ; fi

data_d="/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE"
combine_src="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/extras/20201026_exome_gwas_parallel/4_combine_output.sh"

############################################################
# tmp dir
############################################################
tmp_dir_root="$LOCAL_SCRATCH"
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
# echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT
############################################################

if [ ! -d ${data_d}/${pop}/logs ] ; then mkdir -p ${data_d}/${pop}/logs ; fi

log_template_f=${data_d}/${pop}-batch/logs/ukb24983_exomeOQFE.batch__BATCH__.${trait}.glm.log
log_combined_f=${data_d}/${pop}/logs/ukb24983_exomeOQFE.${trait}.glm.log.gz
plink_template_f=${data_d}/${pop}-batch/ukb24983_exomeOQFE.batch__BATCH__.${trait}.glm.logistic.hybrid
plink_combined_f=${data_d}/${pop}/ukb24983_exomeOQFE.${trait}.glm.logistic.hybrid.gz

log_combined_tmp_f=${tmp_dir}/$(basename ${log_combined_f})
plink_combined_tmp_f=${tmp_dir}/$(basename ${plink_combined_f})

# plink file
if [ ! -s ${plink_combined_f} ] ; then
    bash ${combine_src} \
    --cores ${cores} \
    --n_batch 100 \
    ${plink_template_f} ${plink_combined_tmp_f}
    cp ${plink_combined_tmp_f} ${plink_combined_f}
    echo ${plink_combined_f}
fi                                

# log file
if [ ! -s ${log_combined_f} ] ; then
    bash ${combine_src} \
    --log \
    --cores ${cores} \
    --n_batch 100 \
    ${log_template_f} ${log_combined_tmp_f}
    cp ${log_combined_tmp_f} ${log_combined_f}
    echo ${log_combined_f}
fi

exit 0
########################################################################
bash 2_combine_output_BIN.sh african BIN19 6

#!/bin/bash
set -beEuo pipefail

logit_zst_f=$1
# logit_zst_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/white_british-batch/ukb24983_exomeOQFE.batch76.BIN_FC10002654.PHENO1.glm.logistic.hybrid.zst

plink_log_f=${logit_zst_f%.PHENO1.glm.logistic.hybrid.zst}.log

if [ ! -s ${plink_log_f} ] ; then
    exit 0
fi

is_finished_flag=$(cat ${plink_log_f} | grep 'End time' | wc -l)

logit_zst_renamed=$(echo ${logit_zst_f} | sed -e "s/.PHENO1//g")
plink_log_renamed=$(dirname ${plink_log_f})/logs/$(basename ${plink_log_f%.log}).glm.log

if [ ${is_finished_flag} -eq 1 ] ; then
    echo mv $logit_zst_f $logit_zst_renamed
    echo mv $plink_log_f $plink_log_renamed
    mv $logit_zst_f $logit_zst_renamed
    mv $plink_log_f $plink_log_renamed
fi

exit 0
##########################################

bash 2b_patch_file_move_zst.sh /scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/white_british-batch/ukb24983_exomeOQFE.batch76.BIN_FC10002654.PHENO1.glm.logistic.hybrid.zst

bash 2b_patch_file_move_zst.sh /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_scg/white_british-batch/ukb24983_exomeOQFE.batch7.cancer1050.glm.logistic.hybrid.zst

find /scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/*-batch -maxdepth 1 -type f -name "*.glm.logistic.hybrid.zst" \
| parallel -j6 --eta 'bash 2b_patch_file_move_zst.sh {}'

find /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_scg/*-batch -maxdepth 1 -type f -name "*.glm.logistic.hybrid.zst" \
| parallel -j6 --eta 'bash 2b_patch_file_move_zst.sh {}'
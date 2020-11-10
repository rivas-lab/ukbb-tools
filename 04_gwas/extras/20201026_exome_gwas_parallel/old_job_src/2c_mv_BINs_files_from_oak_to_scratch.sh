#!/bin/bash
set -beEuo pipefail

log_f=$1
# log_f=/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_scg/white_british-batch/logs/ukb24983_exomeOQFE.batch35.cancer1018.glm.log
plink=$(dirname $(dirname ${log_f}))/$(basename ${log_f%.log}).logistic.hybrid
log_f_scratch=$(echo ${log_f} | sed -e "s%/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_scg%/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE%g")
plink_scratch=$(echo ${plink} | sed -e "s%/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_scg%/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE%g")

cp_if_dst_does_not_exist () {
    src=$1
    dst=$2
    if [ ! -f ${dst} ] ; then 
        echo cp ${src} ${dst}
        cp ${src} ${dst}
    fi
}

if [ -f ${plink}.zst ] ; then
    cp_if_dst_does_not_exist ${plink}.zst ${plink_scratch}.zst
else
    cp_if_dst_does_not_exist ${plink} ${plink_scratch}
fi

cp_if_dst_does_not_exist ${log_f} ${log_f_scratch}

exit 0
####################

bash 2c_mv_BINs_files_from_oak_to_scratch.sh /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_scg/white_british-batch/logs/ukb24983_exomeOQFE.batch35.cancer1018.glm.log

find /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_scg/white_british-batch/logs/ -name "*.glm.log" | grep -v QT_ALL \
| parallel --eta -j6 -k 'bash 2c_mv_BINs_files_from_oak_to_scratch.sh {}'


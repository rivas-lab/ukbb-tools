#!/bin/bash
set -beEuo pipefail

source 0_functions.sh

log_f=$1

GBE_ID=$(basename ${log_f} | awk -v FS='.' '{print $2}')

plink_f=$(dirname $(dirname ${log_f}))/$(basename ${log_f%.glm.log.gz}).$(get_plink_suffix ${GBE_ID}).gz

if [ ! -s $plink_f ] ; then
    echo $log_f
    rm $log_f
fi

exit 0
############################

bash 5e_remove_merged_log_if_plink_does_not_exist.sh /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/african/logs/ukb24983_exomeOQFE.HC314.glm.log.gz

find /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/*/logs/ -name "*.log.gz" \
| egrep -v QT_ALL | sort -V | while read f ; do bash 5e_remove_merged_log_if_plink_does_not_exist.sh $f ; done > 5e_remove_merged_log_if_plink_does_not_exist.$(date +%Y%m%d-%H%M%S).txt

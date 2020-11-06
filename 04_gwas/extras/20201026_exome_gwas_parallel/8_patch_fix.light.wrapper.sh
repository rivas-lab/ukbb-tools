#!/bin/bash
set -beEuo pipefail

in_f=$1
shift 1

log_f=/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/extras/20201026_exome_gwas_parallel/logs/patch_fix/$(basename $(dirname ${in_f}))/$(basename ${in_f%.gz}).log

if [ ! -d $(dirname ${log_f}) ] ; then mkdir -p $(dirname ${log_f}) ; fi

echo '---------------' >> ${log_f}
bash 8_patch_fix.light.sh $in_f $@ 2>&1 | tee -a ${log_f}

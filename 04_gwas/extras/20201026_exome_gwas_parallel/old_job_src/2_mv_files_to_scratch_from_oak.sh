#!/bin/bash
set -beEuo pipefail

batch_idx=$1

find "/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_scg/white_british-batch" -name "ukb24983_exomeOQFE.batch${batch_idx}.*.glm.linear*" -type f |
while read f ; do mv $f /scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/white_british-batch/$(basename $f) ; done

exit 0
#############
This script was used to move summary statistics file (WB QT_ALL phenotypes) from OAK (the results from scg4) to scratch

seq 251 500 | parallel -j6 --eta 'bash 2_copy_files_to_scratch.sh {}'


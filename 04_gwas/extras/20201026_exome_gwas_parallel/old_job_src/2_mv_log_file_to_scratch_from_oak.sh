#!/bin/bash
set -beEuo pipefail

find "/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_scg/white_british-batch/logs" -name "ukb24983_exomeOQFE.*QT_ALL.glm.log" -type f |
while read f ; do mv $f /scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/white_british-batch/logs/$(basename $f) ; done

exit 0
##############
This script was used to move log file from QT (QT_ALL) run for WB



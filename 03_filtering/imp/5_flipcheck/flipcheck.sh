#!/bin/bash
set -beEuo pipefail

bash \
~/repos/rivas-lab/ukbb-tools/04_gwas/flipfix/flipcheck.sh \
/oak/stanford/groups/mrivas/ukbb24983/array_imp_combined/pgen/ukb24983_ukb24983_cal_hla_cnv_imp.pvar \
| awk '(NR == 1 || toupper($NF) != toupper($4))' \
| awk '$4 != "P" && $4 != "N"'

#!/bin/bash
set -beEuo pipefail

dir="/oak/stanford/groups/mrivas/ukbb24983/array_imp_combined/pgen_v2"
exts=(
bed
bim
fam
hh
nosex
pgen
pgen.log
psam
pvar.zst
shortnames.bim
)

for ext in ${exts[@]} ; do
    mv ${dir}/ukb24983_ukb24983_cal_hla_cnv_imp.${ext} ${dir}/ukb24983_cal_hla_cnv_imp.${ext}
done
#!/bin/bash

root_dir="/oak/stanford/groups/mrivas/ukbb24983/cal/pgen"

ml load plink

plink --bfile /oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_cALL_v2 --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe --remove /oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_v2.not_used_in_pca.phe --freq

awk '{if ($1==23) {$1="X"; print} else {print}}' plink.frq |  awk '{if ($1==24) {$1="Y"; print} else {print}}' | awk '{if ($1==25) {$1="XY"; print} else {print}}' | awk '{if ($1==26) {$1="MT"; print} else {print}}' | awk '{if (($1=="X") || ($1=="Y") || ($1=="XY") || ($1=="MT")) { print $1"\t"$3"\t"$4"\t"$5}}' >cal_chr_ref_alt_maf.bim

# Grab POS, REF, ALT from .bim
paste <(cut -f1 cal_chr_ref_alt_maf.bim) <(cat /oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_chrX_v2.pvar /oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_chrY_v2.pvar /oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_chrXY_v2.pvar /oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_chrMT_v2.pvar | cut -f4-6) <(cut -f4 cal_chr_ref_alt_maf.bim) | sed '1iCHROM\tPOS\tREF\tALT\tMAF'  > tmp && mv tmp ukb_cal_x_y_xy_mt.bim

rm plink*
rm cal_chr_ref_alt_maf.bim

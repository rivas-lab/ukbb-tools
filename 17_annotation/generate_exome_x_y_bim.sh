#!/bin/bash

ml load plink

plink --bfile /oak/stanford/groups/mrivas/ukbb24983/exome/pgen/spb/data/ukb_exm_spb --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe --remove /oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_v2.not_used_in_pca.phe --freq

awk '{if ($1==23) {$1="X"; print} else {print}}' plink.frq |  awk '{if ($1==24) {$1="Y"; print} else {print}}' | awk '{if (($1=="X") || ($1=="Y")) { print $1"\t"$3"\t"$4"\t"$5}}' >exm_spb_chr_ref_alt_maf.bim

# Grab POS
paste <(cut -f1 exm_spb_chr_ref_alt_maf.bim) <(awk -F'\t' '{if (($1=="X") || ($1=="Y")) {print}}' /oak/stanford/groups/mrivas/ukbb24983/exome/pgen/spb/data/ukb_exm_p_spb.pvar | cut -f2,4,5) <(cut -f4 exm_spb_chr_ref_alt_maf.bim) | sed '1iCHROM\tPOS\tREF\tALT\tMAF'  > tmp && mv tmp ukb_exm_spb_x_y.bim

rm plink*
rm exm_spb_chr_ref_alt_maf.bim

#!/bin/bash

#Concatenate all individuals, put filename in front, get 
find /oak/stanford/groups/mrivas/ukbb24983/cal/mrp/ -name '*_gene.tsv' | xargs -I {} awk -F'\t' 'BEGIN{OFS="\t"}{print FILENAME,$0}' {} | grep -v log_10_BF_sigma_m_var_sim_proximal_coding | cut -d/ -f12 >array_list.tsv

paste <(cut -f1 array_list.tsv | sed 's/_gene.tsv//g') <(cut -f2- array_list.tsv)>array_tmp_list.tsv

sort -k1,1 array_tmp_list.tsv >sorted_array_list.tsv

join -1 1 -2 1 <(cut -f1-2 ../05_gbe/phenotype_info.tsv | sort -k1,1) <(cat sorted_array_list.tsv) >array_mrp_list.tsv

cat <(head -1 /oak/stanford/groups/mrivas/ukbb24983/cal/mrp/extras/adjusted_biomarkers/white_british/BIN10030500_gene.tsv) array_mrp_list.tsv >/oak/stanford/groups/mrivas/ukbb24983/cal/mrp/cal_mrp.tsv

find /oak/stanford/groups/mrivas/ukbb24983/exome/mrp/ -name '*_gene.tsv' | xargs -I {} awk -F'\t' 'BEGIN{OFS="\t"}{print FILENAME,$0}' {} | grep -v log_10_BF_sigma_m_var_sim_proximal_coding | cut -d/ -f12 >exome_list.tsv

paste <(cut -f1 exome_list.tsv | sed 's/_gene.tsv//g') <(cut -f2- exome_list.tsv)>exome_tmp_list.tsv

sort -k1,1 exome_tmp_list.tsv >sorted_exome_list.tsv

join -1 1 -2 1 <(cut -f1-2 ../05_gbe/phenotype_info.tsv | sort -k1,1) <(cat sorted_exome_list.tsv) >exome_mrp_list.tsv

cat <(head -1 /oak/stanford/groups/mrivas/ukbb24983/exome/mrp/extras/adjusted_biomarkers/white_british/BIN10030500_gene.tsv) exome_mrp_list.tsv >/oak/stanford/groups/mrivas/ukbb24983/exome/mrp/exome_mrp.tsv

rm *_list.tsv

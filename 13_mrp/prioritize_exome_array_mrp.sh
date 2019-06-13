find /oak/stanford/groups/mrivas/ukbb24983/cal/mrp/ -name '*_gene.tsv' | xargs -I {} awk -F'\t' 'BEGIN{OFS="\t"}{print FILENAME,$0}' {} | awk -F'\t' '{if ((NF == 22) && ($3 >= 3 || $4 >= 3 || $9 >= 3 || $10 >= 3 || $15 >= 3 || $16 >= 3)) {print $1"\t"$2"\t"$3"\t"$4"\t"$9"\t"$10"\t"$15"\t"$16;} else if ((NF == 16) && ($3 >= 3 || $4 >= 3 || $9 >= 3 || $10 >= 3)) {print $1"\t"$2"\t"$3"\t"$4"\t"$9"\t"$10"\t"NA"\t"NA;} else if ((NF == 10) && ($3 >= 3 || $4 >= 3)) {print $1"\t"$2"\t"$3"\t"$4"\t"NA"\t"NA"\t"NA"\t"NA;}}' | grep -v log_10_BF_sigma_m_var_sim_proximal_coding | cut -d/ -f12 >array_priority_list.txt

paste <(cut -f1 array_priority_list.txt | sed 's/_gene.tsv//g') <(cut -f2- array_priority_list.txt)>array_priority_list.tsv

sort -k1,1 array_priority_list.tsv >sorted_array_priority_list.tsv

join -1 1 -2 1 <(cut -f1-2 ../05_gbe/phenotype_info.tsv | sort -k1,1) <(cat sorted_array_priority_list.tsv) >array_priority_list.tsv

find /oak/stanford/groups/mrivas/ukbb24983/exome/mrp/ -name '*_gene.tsv' | xargs -I {} awk -F'\t' 'BEGIN{OFS="\t"}{print FILENAME,$0}' {} | awk -F'\t' '{if ((NF == 22) && ($3 >= 3 || $4 >= 3 || $9 >= 3 || $10 >= 3 || $15 >= 3 || $16 >= 3)) {print $1"\t"$2"\t"$3"\t"$4"\t"$9"\t"$10"\t"$15"\t"$16;} else if ((NF == 16) && ($3 >= 3 || $4 >= 3 || $9 >= 3 || $10 >= 3)) {print $1"\t"$2"\t"$3"\t"$4"\t"$9"\t"$10"\t"NA"\t"NA;} else if ((NF == 10) && ($3 >= 3 || $4 >= 3)) {print $1"\t"$2"\t"$3"\t"$4"\t"NA"\t"NA"\t"NA"\t"NA;}}' | grep -v log_10_BF_sigma_m_var_sim_proximal_coding | cut -d/ -f12 >exome_priority_list.txt

paste <(cut -f1 exome_priority_list.txt | sed 's/_gene.tsv//g') <(cut -f2- exome_priority_list.txt)>exome_priority_list.tsv

sort -k1,1 exome_priority_list.tsv >sorted_exome_priority_list.tsv

join -1 1 -2 1 <(cut -f1-2 ../05_gbe/phenotype_info.tsv | sort -k1,1) <(cat sorted_exome_priority_list.tsv) >exome_priority_list.tsv

rm array_priority_list.txt exome_priority_list.txt sorted_array_priority_list.tsv sorted_exome_priority_list.tsv
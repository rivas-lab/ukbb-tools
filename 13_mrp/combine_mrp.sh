#!/bin/bash

# sbatch -p mrivas --mem=64000 combine_mrp.sh

dir="/oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_rv_ma_array"

for file in $(ls $dir | grep 0.01 | head -1); do
    zcat $dir/$file | awk -F'\t' '{if (NR == 1) {print "#GBE_ID\tGBE_short_name\tpops\tnum_pops\t"$0}}' > mrp_rv_ma_array_multipop_gbe.tsv
done

for phen in $(cut -f2 ../05_gbe/exome/200k/icdinfo.shortnames.exome.tsv | tail -n +2 | sort -u); do
    if ! grep -Fxq "$phen" gbe_blacklist.tsv; then
        shortname=$(awk -F'\t' -v gbe="$phen" '{if ($2 == gbe) {print $5}}' ../05_gbe/exome/200k/icdinfo.shortnames.exome.tsv)
        file=$(find $dir -name "*_${phen}_*" | grep 0.01)
        numpops=0
        if [ -f "$file" ]; then
            pops=$(echo $(basename $file) | awk -F"_$phen" '{print $1}')
            for pop in white_british non_british_white e_asian s_asian african related others; do
                if [[ "$pops" == *"$pop"* ]]; then
                    numpops=$((numpops+1))
                fi
            done
            if [ $numpops -gt 1 ]; then
                zcat $file | awk -F'\t' -v filename="$phen" -v pops="$pops" -v numpops="$numpops" -v shortname="$shortname" '{if (NR != 1) {print filename"\t"shortname"\t"pops"\t"numpops"\t"$0}}' >> mrp_rv_ma_array_multipop_gbe.tsv
            fi
        fi
    fi
done

zcat /oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_rv_ma_array/white_british_BIN1020483_gene_maf_0.01_se_0.2.tsv.gz | awk -F'\t' '{if (NR == 1) {print "#GBE_ID\tGBE_short_name\tpops\tnum_pops\t"$0}}' > mrp_rv_ma_array_singlepop_gbe.tsv

for phen in $(cut -f1 ../05_gbe/exome/200k/exome_phenotype_info.tsv | tail -n +2 | sort -u); do
    if ! grep -Fxq "$phen" gbe_blacklist.tsv; then
        shortname=$(awk -F'\t' -v gbe="$phen" '{if ($2 == gbe) {print $5}}' ../05_gbe/exome/200k/icdinfo.shortnames.exome.tsv)
        file=$(find $dir -name "*_${phen}_*" | grep 0.01)
        numpops=0
        if [ -f "$file" ]; then
            pops=$(echo $(basename $file) | awk -F"_$phen" '{print $1}')
            for pop in white_british non_british_white e_asian s_asian african related others; do
                if [[ "$pops" == *"$pop"* ]]; then
                    numpops=$((numpops+1))
                fi
            done
            if [ $numpops -eq 1 ]; then
                zcat $file | awk -F'\t' -v filename="$phen" -v pops="$pops" -v numpops="$numpops" -v shortname="$shortname" '{if (NR != 1) {print filename"\t"shortname"\t"pops"\t"numpops"\t"$0}}' >> mrp_rv_ma_array_singlepop_gbe.tsv
            fi
        fi
    fi
done

dir="/oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_rv_exome"

zcat /oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_rv_exome/white_british_INI10030500_gene_maf_0.01_se_100.0.tsv.gz | awk -F'\t' '{if (NR == 1) {print "#GBE_ID\tGBE_short_name\tpops\tnum_pops\t"$0}}' > biomarkers_mrp_rv_exome_gbe.tsv
for phen in $(grep adjusted ../05_gbe/exome/200k/exome_phenotype_info.tsv | cut -f1 | grep -v BIN | sort -u); do
    if ! grep -Fxq "$phen" gbe_blacklist.tsv; then
        shortname=$(awk -F'\t' -v gbe="$phen" '{if ($2 == gbe) {print $5}}' ../05_gbe/exome/200k/icdinfo.shortnames.exome.tsv)
        file=$(find $dir -name "*_${phen}_*" | grep 0.01)
        pops="white_british"
        numpops=1
        zcat $file | awk -F'\t' -v filename="$phen" -v pops="$pops" -v numpops="$numpops" -v shortname="$shortname" '{if (NR != 1) {print filename"\t"shortname"\t"pops"\t"numpops"\t"$0}}' >> biomarkers_mrp_rv_exome_gbe.tsv
    fi
done

dir="/oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_rv_exome_var"

zcat /oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_rv_exome_var/white_british_INI10030500_gene_maf_0.01_se_100.0.tsv.gz | awk -F'\t' '{if (NR == 1) {print "#GBE_ID\tGBE_short_name\tpops\tnum_pops\t"$0}}' > biomarkers_mrp_rv_exome_var_gbe.tsv
for phen in $(grep adjusted ../05_gbe/exome/200k/exome_phenotype_info.tsv | cut -f1 | grep -v BIN | sort -u); do
    if ! grep -Fxq "$phen" gbe_blacklist.tsv; then
        shortname=$(awk -F'\t' -v gbe="$phen" '{if ($2 == gbe) {print $5}}' ../05_gbe/exome/200k/icdinfo.shortnames.exome.tsv)
        file=$(find $dir -name "*_${phen}_*" | grep 0.01)
        pops="white_british"
        numpops=1
        zcat $file | awk -F'\t' -v filename="$phen" -v pops="$pops" -v numpops="$numpops" -v shortname="$shortname" '{if (NR != 1) {print filename"\t"shortname"\t"pops"\t"numpops"\t"$0}}' >> biomarkers_mrp_rv_exome_var_gbe.tsv
    fi
done


#for file in $(ls $dir | head -1); do
#    awk -F'\t' -v filename="$name" '{if (NR == 1) {print "GBE_ID\t"$0}}' $dir/$file > mrp_rv_array_gbe.tsv
#done


dir="/oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_rv_ma_exome"

for file in $(ls $dir | grep 0.01 | head -1); do
    zcat $dir/$file | awk -F'\t' '{if (NR == 1) {print "#GBE_ID\tGBE_short_name\tpops\tnum_pops\t"$0}}' > mrp_rv_ma_exome_multipop_gbe.tsv
done

for phen in $(cut -f2 ../05_gbe/exome/200k/icdinfo.shortnames.exome.tsv | tail -n +2 | sort -u); do
    if ! grep -Fxq "$phen" gbe_blacklist.tsv; then
        shortname=$(awk -F'\t' -v gbe="$phen" '{if ($2 == gbe) {print $5}}' ../05_gbe/exome/200k/icdinfo.shortnames.exome.tsv)
        file=$(find $dir -name "*_${phen}_*" | grep 0.01)
        numpops=0
        if [ -f "$file" ]; then
            pops=$(echo $(basename $file) | awk -F"_$phen" '{print $1}')
            for pop in white_british non_british_white e_asian s_asian african related others; do
                if [[ "$pops" == *"$pop"* ]]; then
                    numpops=$((numpops+1))
                fi
            done
            if [ $numpops -gt 1 ]; then
                zcat $file | awk -F'\t' -v filename="$phen" -v pops="$pops" -v numpops="$numpops" -v shortname="$shortname" '{if (NR != 1) {print filename"\t"shortname"\t"pops"\t"numpops"\t"$0}}' >> mrp_rv_ma_exome_multipop_gbe.tsv
            fi
        fi
    fi
done

zcat /oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_rv_ma_exome/white_british_BIN_FC20004825_gene_maf_0.01_se_0.2.tsv.gz | awk -F'\t' '{if (NR == 1) {print "#GBE_ID\tGBE_short_name\tpops\tnum_pops\t"$0}}' > mrp_rv_ma_exome_singlepop_gbe.tsv

for phen in $(cut -f1 ../05_gbe/exome/200k/exome_phenotype_info.tsv | tail -n +2 | sort -u); do
    if ! grep -Fxq "$phen" gbe_blacklist.tsv; then
        shortname=$(awk -F'\t' -v gbe="$phen" '{if ($2 == gbe) {print $5}}' ../05_gbe/exome/200k/icdinfo.shortnames.exome.tsv)
        file=$(find $dir -name "*_${phen}_*" | grep 0.01)
        numpops=0
        if [ -f "$file" ]; then
            pops=$(echo $(basename $file) | awk -F"_$phen" '{print $1}')
            for pop in white_british non_british_white e_asian s_asian african related others; do
                if [[ "$pops" == *"$pop"* ]]; then
                    numpops=$((numpops+1))
                fi
            done
            if [ $numpops -eq 1 ]; then
                zcat $file | awk -F'\t' -v filename="$phen" -v pops="$pops" -v numpops="$numpops" -v shortname="$shortname" '{if (NR != 1) {print filename"\t"shortname"\t"pops"\t"numpops"\t"$0}}' >> mrp_rv_ma_exome_singlepop_gbe.tsv
            fi
        fi
    fi
done

python combine_mrp.py
rm mrp_rv_ma_exome_singlepop_gbe.tsv
rm mrp_rv_ma_exome_multipop_gbe.tsv
rm mrp_rv_ma_array_singlepop_gbe.tsv
rm mrp_rv_ma_array_multipop_gbe.tsv
bgzip -f mrp_rv_ma_exome_gbe.tsv
bgzip -f mrp_rv_ma_array_gbe.tsv
bgzip -f biomarkers_mrp_rv_ma_exome_gbe.tsv
bgzip -f biomarkers_mrp_rv_ma_array_gbe.tsv
bgzip -f biomarkers_mrp_rv_exome_gbe.tsv
bgzip -f biomarkers_mrp_rv_exome_var_gbe.tsv

#!/bin/bash

#dir="/oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_rv_array"

#for file in $(ls $dir | head -1); do
#    awk -F'\t' -v filename="$name" '{if (NR == 1) {print "GBE_ID\t"$0}}' $dir/$file > mrp_rv_array_gbe.txt
#done

#for phen in $(cut -f1 ../05_gbe/phenotype_info.tsv | tail -n +2); do
#    file=$(find $dir -name "*_${phen}_*")
#    if [ -f "$file" ]; then
#        awk -F'\t' -v filename="$phen" '{if (NR != 1) {print filename"\t"$0}}' $file >> mrp_rv_array_gbe.txt
#    fi
#done

dir="/oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_rv_ma_array"

for file in $(ls $dir | head -1); do
    awk -F'\t' '{if (NR == 1) {print "#GBE_ID\tpops\t"$0}}' $dir/$file > mrp_rv_ma_array_multipop_gbe.txt
done

for phen in $(cut -f1 ../05_gbe/array-combined/phenotype_info.tsv | tail -n +2); do
#for phen in INI22154; do
    file=$(find $dir -name "*_${phen}_*")
    numpops=0
    if [ -f "$file" ]; then
        pops=$(echo $(basename $file) | awk -F"_$phen" '{print $1}')
        for pop in white_british non_british_white e_asian s_asian african related others; do
            if [[ "$pops" == *"$pop"* ]]; then
                numpops=$((numpops+1))
            fi
        done
        echo $numpops
        echo $file
        if [ $numpops -gt 1 ]; then
            awk -F'\t' -v filename="$phen" -v pops="$pops" '{if (NR != 1) {print filename"\t"pops"\t"$0}}' $file >> mrp_rv_ma_array_multipop_gbe.txt
        fi
    fi
done

awk -F'\t' '{if (NR == 1) {print "#GBE_ID\tpops\t"$0}}' /oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_rv_ma_array/white_british_BIN20474_gene_maf_0.01_se_0.2.tsv > mrp_rv_ma_array_singlepop_gbe.txt

for phen in $(cut -f1 ../05_gbe/array-combined/phenotype_info.tsv | tail -n +2); do
#for phen in INI22154; do
    file=$(find $dir -name "*_${phen}_*")
    numpops=0
    if [ -f "$file" ]; then
        pops=$(echo $(basename $file) | awk -F"_$phen" '{print $1}')
        for pop in white_british non_british_white e_asian s_asian african related others; do
            if [[ "$pops" == *"$pop"* ]]; then
                numpops=$((numpops+1))
            fi
        done
        echo $numpops
        echo $file
        if [ $numpops -eq 1 ]; then
            awk -F'\t' -v filename="$phen" -v pops="$pops" '{if (NR != 1) {print filename"\t"pops"\t"$0}}' $file >> mrp_rv_ma_array_singlepop_gbe.txt
        fi
    fi
done

python combine_mrp.py
rm mrp_rv_ma_array_singlepop_gbe.txt 
rm mrp_rv_ma_array_multipop_gbe.txt 

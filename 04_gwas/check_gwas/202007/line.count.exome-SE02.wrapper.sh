#!/bin/bash

for pop in african e_asian s_asian non_british_white white_british related others; do
#for pop in metal; do
    echo $pop
    ls /scratch/groups/mrivas/ukbb24983/exome/gwas-qc-SE02/${pop}/*tsv | cut -d'.' -f2 >complete
    for gbe_id in $(grep -Fvwf complete <(ls /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/current/${pop}/ | grep -v logs | cut -d'.' -f2 | sort -u)); do
        echo $gbe_id
        exome_f=$(find /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/current/${pop}/ -name "*.$gbe_id.*" | grep -v exome-spb.log)
        #exome_f=/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/current/${pop}/${gbe_id}.metal.tsv.gz
        echo $exome_f
        if [ -f $exome_f ]; then
            sbatch -p owners,mrivas line.count.exome-SE02.sh $exome_f
        else
            echo "$gbe_id not found."
        fi
    done
done

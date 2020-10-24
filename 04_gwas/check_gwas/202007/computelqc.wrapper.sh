#!/bin/bash

#for pop in african e_asian s_asian white_british non_british_white related others metal; do
for pop in metal; do
    echo $pop
    ls /scratch/groups/mrivas/ukbb24983/array-combined/gwas-qc/${pop}/*txt | cut -d'.' -f2 >complete
    for gbe_id in $(grep -Fwvf complete <(ls /oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current/${pop}/ | grep -v logs | cut -d'.' -f1)); do
        echo $gbe_id
        array_f=/oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current/${pop}/${gbe_id}.metal.tsv.gz
        ls $array_f
        if [ -f $array_f ]; then
            sbatch -p owners,mrivas --mem=16000 line.count.sh $array_f $pop
        else
            echo "$gbe_id not found."
        fi
    done
done

#!/bin/bash

for pop in african e_asian s_asian non_british_white white_british related others; do
#for pop in metal; do
    echo $pop
    ls /scratch/groups/mrivas/ukbb24983/array-combined/gwas-qc-SE02/${pop}/*tsv | cut -d'.' -f2 >complete
    for gbe_id in $(grep -Fvwf complete <(ls /oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current/${pop}/ | grep -v logs | cut -d'.' -f2)); do
    #for gbe_id in $(cat to_update); do
        echo $gbe_id
        array_f=$(find /oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current/${pop}/ -name "*.$gbe_id.*" | grep -v array-combined.log)
        #array_f=/oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current/${pop}/${gbe_id}.metal.tsv.gz
        echo $array_f
        if [ -f $array_f ]; then
            sbatch -p owners,mrivas line.count-SE02.sh $array_f
        else
            echo "$gbe_id not found."
        fi
    done
done

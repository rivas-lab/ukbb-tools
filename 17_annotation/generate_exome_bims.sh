#!/bin/bash


for pop in african e_asian non_british_white s_asian white_british; do
    echo $pop
    zcat /oak/stanford/groups/mrivas/ukbb24983/exome/pgen/spb/data/ukb_exm_spb-${pop}-variant_annots_gbe.tsv.gz | cut -f1-4,18 | tail -n +2 | sed '1iCHROM\tPOS\tREF\tALT\tMAF' > ukb_exm_spb-${pop}.bim
done

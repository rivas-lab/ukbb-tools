#!/bin/bash

touch combined_sites.txt

for file in /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/qt/BBJ.AG.combined.txt.gz /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/bin/BBJ.AF.txt.gz /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/bin/BBJ.colorectal_cancer.txt.gz /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/bin/BBJ.POAG.autosome.txt.gz /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/bin/BBJ.POAG.chrX.female.txt.gz /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/bin/BBJ.POAG.chrX.male.txt.gz /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/bin/BBJ.RA.transethnic.combined.txt.gz /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/bin/BBJ.smoking_stop.autosome.txt.gz /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/bin/BBJ.smoking_stop.combined.female.txt.gz /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/bin/BBJ.smoking_stop.combined.male.txt.gz /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/bin/BBJ.T2D.txt.gz; do
    cat combined_sites.txt <(zcat $file | cut -f1,2,4,5) | sort -u >combined_sites_tmp.txt && mv combined_sites_tmp.txt combined_sites.txt
done

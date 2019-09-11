#!/bin/bash

touch bbj_combined_sites.txt

for file in /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/bin/BBJ.AF.txt.gz /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/bin/BBJ.colorectal_cancer.txt.gz /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/bin/BBJ.POAG.autosome.txt.gz /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/bin/BBJ.POAG.chrX.female.txt.gz /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/bin/BBJ.POAG.chrX.male.txt.gz /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/bin/BBJ.RA.transethnic.combined.txt.gz /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/bin/BBJ.smoking_stop.autosome.txt.gz /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/bin/BBJ.smoking_stop.combined.female.txt.gz /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/bin/BBJ.smoking_stop.combined.male.txt.gz /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/bin/BBJ.T2D.txt.gz; do
    #echo $file
    #zcat $file | cut -f9 | grep ADD
    cat bbj_combined_sites.txt <(zcat $file | cut -f1,2,4,5,9) | sort -u >bbj_combined_sites_tmp.txt && mv bbj_combined_sites_tmp.txt bbj_combined_sites.txt
done

for file in /oak/stanford/groups/mrivas/bbj/flipfixed_combined_aut_X/qt/BBJ.AG.combined.txt.gz; do
    #echo $file
    #zcat $file | cut -f8 | grep ADD
    cat bbj_combined_sites.txt <(zcat $file | cut -f1,2,4,5,8) | sort -u >bbj_combined_sites_tmp.txt && mv bbj_combined_sites_tmp.txt bbj_combined_sites.txt
done

grep -v e bbj_combined_sites.txt | grep -v CHROM >bbj_combined_sites_tmp.txt && mv bbj_combined_sites_tmp.txt bbj_combined_sites.txt

awk -F'\t' '{if ($5 >= 0.5) { print $1"\t"$2"\t"$3"\t"$4"\t"(1-$5) } else {print $0}}' bbj_combined_sites.txt >bbj_combined_sites_tmp.txt && mv bbj_combined_sites_tmp.txt bbj_combined_sites.txt

sed -i '1s/^/CHROM\tPOS\tREF\tALT\tMAF\n/' bbj_combined_sites.txt

sbatch --mem=64000 -p mrivas --wrap="python compute_maf.py"

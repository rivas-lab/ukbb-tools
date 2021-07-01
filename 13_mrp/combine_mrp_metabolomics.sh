#!/bin/bash

# sbatch -p mrivas --mem=64000 combine_mrp.sh

dir="/oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_metabolomics_mpc_pli_new"

for file in $(ls $dir | grep "0\.01" | head -1); do
    zcat $dir/$file | awk -F'\t' '{if (NR == 1) {print "#GBE_ID\tGBE_short_name\tpops\tnum_pops\t"$0}}' > mrp_rv_metabolomics.tsv
done

for phen in $(cat metabolomics_phenos); do
    if ! grep -Fxq "$phen" gbe_blacklist.tsv; then
        shortname=$(awk -F'\t' -v gbe="$phen" '{if ($2 == gbe) {print $5}}' ../05_gbe/exome/200k/icdinfo.shortnames.exome.tsv)
        file=$(find $dir -name "*_${phen}_*" | grep "0\.01")
        numpops=0
        if [ -f "$file" ]; then
            pops=$(echo $(basename $file) | awk -F"_$phen" '{print $1}')
            for pop in white_british non_british_white e_asian s_asian african related others; do
                if [[ "$pops" == *"$pop"* ]]; then
                    numpops=$((numpops+1))
                fi
            done
            if [ $numpops -gt 1 ]; then
                zcat $file | awk -F'\t' -v filename="$phen" -v pops="$pops" -v numpops="$numpops" -v shortname="$shortname" '{if (NR != 1) {print filename"\t"shortname"\t"pops"\t"numpops"\t"$0}}' >> mrp_rv_metabolomics.tsv
            fi
        fi
    fi
done

python combine_mrp_metabolomics.py
bgzip -f mrp_rv_metabolomics.tsv

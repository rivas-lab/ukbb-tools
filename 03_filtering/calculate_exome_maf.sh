#/bin/bash

ml load plink2

qcDir="/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification"
dataPrefix="/oak/stanford/groups/mrivas/private_data/ukbb/24983/exome/pgen/spb/data/ukb_exm_spb"

for pop in white_british non_british_white african s_asian e_asian; do
    plink2 --keep $qcDir/ukb24983_${pop}.phe --freq;
done

#!/bin/bash
set -beEuo pipefail

# for pop in 'african' 's_asian' 'e_asian' 'related' 'others' 'non_british_white' 'white_british' ; do
for pop in 'african' 's_asian' 'e_asian' 'related' 'others' 'non_british_white' ; do
#for pop in 'white_british' ; do
jn=gwas.QTs.${pop}

sbatch $([ "${pop}" == "white_british" ] && echo "-p mrivas --qos=high_p" || echo "-p normal,mrivas") \
    --nodes=1 --mem=7500 --cores=2 --time=12:00:00 \
    --job-name=${jn} --output=logs/${jn}.%A_%a.out --error=logs/${jn}.%A_%a.err \
    --array=$([ "${pop}" == "white_british" ] && echo "2-1000" || echo "1-100") \
    1_plink.gwas.QT_ALL.sh ${pop}
done


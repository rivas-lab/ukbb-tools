#!/bin/bash
set -beEuo pipefail

GBE_ID=$1

# for pop in 'african' 's_asian' 'e_asian' 'related' 'others' 'non_british_white' 'white_british' ; do
# for pop in 'african' 's_asian' 'e_asian' 'related' 'others' 'non_british_white' ; do
for pop in 'white_british' ; do
jn=gwas.${GBE_ID}.${pop}

sbatch $([ "${pop}" == "white_british" ] && echo "-p mrivas --qos=high_p" || echo "-p owners") \
    --nodes=1 --mem=9000 --cores=2 --time=12:00:00 \
    --job-name=${jn} --output=logs_scratch/${jn}.%A_%a.out --error=logs_scratch/${jn}.%A_%a.err \
    --array=$([ "${pop}" == "white_british" ] && echo "1-100" || echo "1-100") \
    1_plink.gwas.GBE_ID.sh ${pop} ${GBE_ID}
done


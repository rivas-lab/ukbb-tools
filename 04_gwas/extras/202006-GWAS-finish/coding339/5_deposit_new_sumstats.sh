#!/bin/bash
set -beEuo pipefail

pop=$1

gwas_dir="/oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas"

GBE_IDs=(
    'INI21049' 'INI21051' 'INI21052'
    'INI21053' 'INI21054' 'INI21055' 'INI21056'
    'INI21058' 'INI21059' 'INI21060' 'INI21061'
)
current_d="${gwas_dir}/current/${pop}"
data_d="${gwas_dir}/2005693/37855/${pop}"


for GBE_ID in ${GBE_IDs[@]} ; do
    new_f=data/${pop}/ukb24983_v2_hg19.${GBE_ID}.array-combined.glm.linear.gz
    cp $new_f ${data_d}/
    ln -s ${data_d}/$(basename ${new_f}) ${current_d}/$(basename ${new_f})
done

exit 0
####################
# instructions 2020/7/4
for pop in 'non_british_white' 'african' 's_asian' 'e_asian' 'related' 'others' ; do bash 5_deposit_new_sumstats.sh $pop ; done
bash 5_deposit_new_sumstats.sh white_british

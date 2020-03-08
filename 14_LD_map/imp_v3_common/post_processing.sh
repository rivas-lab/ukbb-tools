#!/bin/bash
set -beEuo pipefail

pops=(
"white_british"
"non_british_white"
"african"
"s_asian"
"e_asian"
)

ldmap_post_processing () {
    local pop=$1
    local chr=$2

    local data_dir="/scratch/groups/mrivas/ukbb24983/imp/ldmap_common"
    local basename="${data_dir}/ukb24983_imp_common_chr${chr}_v3.${pop}.ld_map"

    if [ ! -f ${basename}.tsv.gz ] ; then
        zcat ${basename}.ld.gz \
        | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7}'  \
        | sed -e "s/^CHR_A/#CHR_A/g" \
        | bgzip -l9 -@6 > ${basename}.tsv.gz
    fi

    if [ ! -f ${basename}.tsv.gz.tbi ] ; then
        tabix -c '#' -s 1 -b 2 -e 5 ${basename}.tsv.gz
    fi
}

for chr in $(seq 1 22) X XY ; do for pop in "${pops[@]}" ; do
    ldmap_post_processing ${pop} ${chr}
done; done

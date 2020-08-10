#!/bin/bash
set -beEuo pipefail

GBE_ID=HC382
n_batch=100

out_dir="/oak/stanford/groups/mrivas/users/ytanigaw/sandbox/20200725-plink2-firth/"

##################

get_plink_suffix () {
    local GBE_ID=$1

    GBE_CAT=$(echo $GBE_ID | sed -e "s/[0-9]//g")

    if [ "${GBE_CAT}" == "QT_FC" ] || [ "${GBE_CAT}" == "INI" ] ; then
        echo glm.linear
    else
        echo glm.logistic.hybrid
    fi
}

##################

suffix=$(get_plink_suffix ${GBE_ID})

for out_d in ${out_dir}/firth-residualize ${out_dir}/firth-default ${out_dir}/cc-residualize ; do
# combine the log files

seq ${n_batch} | tr ' ' '\n' | while read batch_idx ; do
    echo "## ${out_d}/ukb24983_v2_hg19.${GBE_ID}.batch${batch_idx}.log"
    cat ${out_d}/ukb24983_v2_hg19.${GBE_ID}.batch${batch_idx}.log
done | bgzip -l9 -@6 > ${out_d}/ukb24983_v2_hg19.${GBE_ID}.log.gz

# combine the summary statistics files

{
    cat ${out_d}/ukb24983_v2_hg19.${GBE_ID}.batch1.${suffix} | egrep '^#'

    seq ${n_batch} | tr ' ' '\n' | while read batch_idx ; do
        cat ${out_d}/ukb24983_v2_hg19.${GBE_ID}.batch${batch_idx}.${suffix} | egrep -v '^#'
    done | sort -k1,1V -k2,2n -k3,3
} | bgzip -l9 -@6 > ${out_d}/ukb24983_v2_hg19.${GBE_ID}.${suffix}.gz

done

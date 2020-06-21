#!/bin/bash
set -beEuo pipefail

in_file=$1
out_file=$(dirname $(dirname ${in_file}))/summary_stats_hg19/$(basename ${in_file} .gz).hg19

if [ ! -d $(dirname ${out_file}) ] ; then mkdir -p $(dirname ${out_file}) ; fi

if [ ! -f ${out_file}.gz.tbi ] ; then

    Rscript 4_FinnGenR3_liftOver.R ${in_file} ${out_file}
    {
        cat ${out_file} | egrep '^#'
        cat ${out_file} | egrep -v '^#' | sort --parallel 4 -k1,1V -k2,2n -k5,5 
    }| bgzip -l9 -@4 > ${out_file}.gz
    tabix -f -c '#' -s 1 -b2 -e2 ${out_file}.gz

    echo ${out_file}.gz

fi


#!/bin/bash
set -beEuo pipefail

sumstats_file=$(readlink -f $1)

echo $sumstats_file

if [ ${sumstats_file%.gz}.gz != ${sumstats_file} ] ; then
    bgzip --compress-level 9 -f ${sumstats_file}
    log_file=$(echo ${sumstats_file} | sed -e 's/.linear$//g' | sed -e 's/.logistic.hybrid$//g' | sed -e 's/.glm$/.log/g' )
    mv $log_file $(dirname $log_file)/logs/
fi


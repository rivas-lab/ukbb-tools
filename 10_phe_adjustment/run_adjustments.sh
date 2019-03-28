#!/bin/bash
set -beEuo pipefail

adjustment_tsv="$(dirname $0)/adjustingbiomarker.tsv"
script="$(dirname $0)/phe_adjustment.py"
phe_dir="/oak/stanford/groups/mrivas/dev-ukbb-tools/phenotypes"

in_dir="${phe_dir}/2001440"
out_dir="${phe_dir}/2001440_adjusted"

if [ ! -d ${out_dir} ] ; then mkdir -p ${out_dir} ; fi

cat $adjustment_tsv | egrep -v '^#' | egrep -v '^Intercept' \
| awk '{print $2, $4}' \
| while read GBE mult ; do
    echo $GBE  $mult
    cat ${in_dir}/${GBE}.phe | python $script --mult $mult > ${out_dir}/${GBE}.adjusted.phe
done


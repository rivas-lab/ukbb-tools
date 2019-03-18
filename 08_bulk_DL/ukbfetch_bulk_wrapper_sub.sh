#!/bin/bash
set -beEuo pipefail

bulk_file=$1
batch_idx=$2
out_dir=$3

n_lines=$(cat ${bulk_file} | wc -l)
n_batches=$( perl -e "print(int((${n_lines} + 999)/1000))" )

get_batch_start () {
    batch_id=$1

    echo $( perl -e "print(int((${batch_id} - 1) * 1000 + 1))")
}

get_batch_end () {
    batch_id=$1

    if [ $batch_id -eq ${n_batches} ] ; then
	echo $n_lines
    else
	echo $( perl -e "print(int(${batch_id} * 1000))")
    fi
}

fetch_batch () {
    file=$1
    idx=$2
    out_d=$3

    cat ${file} \
	| awk -v rs=$(get_batch_start $idx) -v re=$(get_batch_end $idx) 'rs <= NR && NR <= re' \
	| ${out_d}/ukbfetch -b/dev/stdin -o$(basename ${file}).bulk${idx}
}

# copy ukbfetch
cd ${out_dir}
if [ ! -f ${out_dir}/ukbfetch ] ; then
	wget  -nd  biobank.ndph.ox.ac.uk/showcase/util/ukbfetch
	chmod 775 ${out_dir}/ukbfetch
fi

# run task
fetch_batch ${bulk_file} ${batch_idx} ${out_dir}


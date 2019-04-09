#!/bin/bash
set -beEuo pipefail

usage () {
    echo "usage: $0 in.bed.lst assembly out.dir"
}

get_file_list () {
    cat $1 | egrep -v '^#' | awk '{print $1}'
}

source $(dirname $(readlink -f $0))/great_misc_func.sh

if [ $# -lt 3 ] ; then usage >&2 ; exit 1 ; fi
in_file_list=$1
assembly=$2
out_d=$3

# create a temp directory
tmp_root="/dev/shm/"
tmp_dir=$(mktemp -d -p ${tmp_root} tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX) ; echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT
for sub_d in great out ; do
    if [ ! -d ${tmp_dir}/${sub_d} ] ; then mkdir -p ${tmp_dir}/${sub_d} ; fi
done


get_file_list ${in_file_list} | while read file ; do
    echo "Running GREAT for ${file} .." >&2
    if [ ! -f $tmp_dir/great/$(basename ${file%.gz} .bed) ] ; then
        bash $(dirname $(readlink -f $0))/great_wrapper.sh $file $assembly $tmp_dir/great/$(basename ${file%.gz} .bed)
    fi
done

for onto in $( get_onto_list $assembly | tr "," " " ) ; do
    echo "Aggregating summary statistics for $onto .." >&2
    {
        echo "#Phe $( great_cols | sed -e 's/#//g' )" | tr " " "\t"
        get_file_list ${in_file_list} | while read file ; do
            cat ${tmp_dir}/great/$(basename ${file%.gz} .bed)/${onto}.tsv | egrep -v '^#' \
            | awk -v file_id=$(basename ${file%.gz} .bed) -v OFS='\t' '{print file_id, $0}'
        done 
    } > ${tmp_dir}/out/${onto}.tsv
    gzip ${tmp_dir}/out/${onto}.tsv
done

cp -ar ${tmp_dir}/out $out_d


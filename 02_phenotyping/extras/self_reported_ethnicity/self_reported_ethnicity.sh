#!/bin/bash
set -beEuo pipefail

in_f="/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/9796/24611/download/ukb24611.tab"
out_f="/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/self_reported_ethnicity/misc/ukb9796_ukb24611_f21000.tsv"

tmp_dir_root=$LOCAL_SCRATCH
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT

tmp_in_f=${tmp_dir}/$(basename $in_f)
tmp_out_f=${tmp_dir}/$(basename $out_f)

date
echo "copy the input file"
cp ${in_f} ${tmp_in_f}

date
echo "extract relevant cols.."
cat ${tmp_in_f} \
    | cut -f 1,3986,3987,3988 \
    | sed -e "s/f\\.//g" \
    | sed -e "s/^eid/IID/g" > ${tmp_out_f}

date
echo "copy the results"
if [ ! -d $(dirname ${out_f}) ] ; then mkdir -p $(dirname ${out_f}) ; fi
cp ${tmp_out_f} ${out_f}

date
echo "done"


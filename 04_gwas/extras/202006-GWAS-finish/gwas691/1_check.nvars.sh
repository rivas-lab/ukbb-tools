#!/bin/bash
# set -beEuo pipefail

wc_l=$1

############################################################
# tmp dir
############################################################
tmp_dir_root="$LOCAL_SCRATCH"
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
# echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT

############################################################

cat ../gwas-current-gz-wc.20200629-131527.tsv | awk -v FS='\t' -v wc_l=${wc_l} '($3==wc_l){print $2}' \
| while read f ; do
    echo $f >&2
    zcat $f | cut -f3 > ${tmp_dir}/$(basename $(dirname $f))_$(basename $f).IDs.txt
done

find ${tmp_dir} -type f -name "*.IDs.txt" | while read f ; do cat $f ; done | sort --parallel 6 -u > 1_check.nvars.out.lst

cat 1_check.nvars.out.lst | wc -l

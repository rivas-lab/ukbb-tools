#!/bin/bash
set -beEuo pipefail

enc_file=$(readlink -f $1)
key_file=$(readlink -f $2)

out_dir=$(dirname ${enc_file})
prefix=$(basename ${enc_file%.enc})

ml load ukbb-showcase-utils R

cd $out_dir

ukbunpack ${prefix}.enc $(cat $key_file | awk 'NR>1' | rev | cut -c2- | rev)

wget -nd biobank.ndph.ox.ac.uk/showcase/util/encoding.ukb

# documentation
ukbconv ${prefix}.enc_ukb docs
mv ${out_dir}/${prefix}.log ${out_dir}/${prefix}.html.log

# R tab file
ukbconv ${prefix}.enc_ukb r
mv ${out_dir}/${prefix}.log ${out_dir}/${prefix}.tab.log

echo ${out_dir}/${prefix}.tab

echo "[$0] tab.columns" >&2
ukbtabcols.sh ${out_dir}/${prefix}.tab > ${out_dir}/${prefix}.tab.columns

echo "[$0] tab.columns.summary.tsv.gz" >&2
bash $(dirname $(readlink -f $0))/compute_col-wise_md5sum.sh ${out_dir}/${prefix}.tab > ${out_dir}/${prefix}.tab.columns.summary.tsv
gzip -9 ${out_dir}/${prefix}.tab.columns.summary.tsv

cd -


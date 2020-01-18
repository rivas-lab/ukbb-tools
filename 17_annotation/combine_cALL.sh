#!/bin/bash
set -beEuo pipefail


data_dir="/oak/stanford/groups/mrivas/users/guhan/repos/ukbb-tools/16_annotation/imp_annots"

{
  ! cat ${data_dir}/ukb24983_imp_chr1_v3.pvar.zst_cf_vep.tsv | head -n1 | sed -e "s/CHROM/#CHROM/g"

for c in $(seq 1 22) ; do
    cat ${data_dir}/ukb24983_imp_chr${c}_v3.pvar.zst_cf_vep.tsv | awk 'NR>1'
done

} | bgzip -@ 4 -l 9 > /oak/stanford/groups/mrivas/ukbb24983/imp/annotation/annot.tsv.gz


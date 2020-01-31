#!/bin/bash

mv /oak/stanford/groups/mrivas/bbj/variant_annotation/bbj_variant_annots.tsv.bgz /oak/stanford/groups/mrivas/bbj/variant_annotation/bbj_variant_annots.tsv.gz

zcat /oak/stanford/groups/mrivas/bbj/variant_annotation/bbj_variant_annots.tsv.gz | cut -f3- >/oak/stanford/groups/mrivas/bbj/variant_annotation/bbj_variant_annots.tsv

#remove all the nulls
sed -i 's/null,//g' /oak/stanford/groups/mrivas/bbj/variant_annotation/bbj_variant_annots.tsv
sed -i 's/,null//g' /oak/stanford/groups/mrivas/bbj/variant_annotation/bbj_variant_annots.tsv
#replace with NAs
sed -i 's/\[null\]/NA/g' /oak/stanford/groups/mrivas/bbj/variant_annotation/bbj_variant_annots.tsv
#get rid of quotes
sed -i 's/"//g' /oak/stanford/groups/mrivas/bbj/variant_annotation/bbj_variant_annots.tsv
#replace NAs with quotes
sed -i 's/\<NA\>/""/g' /oak/stanford/groups/mrivas/bbj/variant_annotation/bbj_variant_annots.tsv
#remove brackets
sed -i 's/\[//g' /oak/stanford/groups/mrivas/bbj/variant_annotation/bbj_variant_annots.tsv
sed -i 's/\]//g' /oak/stanford/groups/mrivas/bbj/variant_annotation/bbj_variant_annots.tsv

awk -F'\t' 'BEGIN{OFS="\t"}{if (NR == 1) { print $0 } else { $13="CSQ="$13; print}}' /oak/stanford/groups/mrivas/bbj/variant_annotation/bbj_variant_annots.tsv >bbj_variant_annots_tmp.tsv && mv bbj_variant_annots_tmp.tsv /oak/stanford/groups/mrivas/bbj/variant_annotation/bbj_variant_annots_tmp.tsv

#gzip -f /oak/stanford/groups/mrivas/bbj/variant_annotation/bbj_variant_annots.tsv

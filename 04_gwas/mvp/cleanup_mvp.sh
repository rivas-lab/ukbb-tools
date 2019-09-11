#!/bin/bash

#mv /oak/stanford/groups/mrivas/public_data/mvp_20190829/mvp_variant_annots.tsv.bgz /oak/stanford/groups/mrivas/public_data/mvp_20190829/mvp_variant_annots.tsv.gz

#zcat /oak/stanford/groups/mrivas/public_data/mvp_20190829/mvp_variant_annots.tsv.gz | cut -f3- >/oak/stanford/groups/mrivas/public_data/mvp_20190829/mvp_variant_annots.tsv

#remove all the nulls
#sed -i 's/null,//g' /oak/stanford/groups/mrivas/public_data/mvp_20190829/mvp_variant_annots.tsv
#sed -i 's/,null//g' /oak/stanford/groups/mrivas/public_data/mvp_20190829/mvp_variant_annots.tsv
#replace with NAs
#sed -i 's/\[null\]/NA/g' /oak/stanford/groups/mrivas/public_data/mvp_20190829/mvp_variant_annots.tsv
#get rid of quotes
#sed -i 's/"//g' /oak/stanford/groups/mrivas/public_data/mvp_20190829/mvp_variant_annots.tsv
#replace NAs with quotes
#sed -i 's/\<NA\>/""/g' /oak/stanford/groups/mrivas/public_data/mvp_20190829/mvp_variant_annots.tsv
#remove brackets
#sed -i 's/\[//g' /oak/stanford/groups/mrivas/public_data/mvp_20190829/mvp_variant_annots.tsv
#sed -i 's/\]//g' /oak/stanford/groups/mrivas/public_data/mvp_20190829/mvp_variant_annots.tsv
awk -F'\t' 'BEGIN{OFS="\t"}{if (NR == 1) { print $0 } else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\tCSQ="$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30}}' /oak/stanford/groups/mrivas/public_data/mvp_20190829/mvp_variant_annots.tsv >mvp_variant_annots_tmp.tsv #&& mv mvp_variant_annots_tmp.tsv /oak/stanford/groups/mrivas/public_data/mvp_20190829/mvp_variant_annots.tsv

#gzip -f /oak/stanford/groups/mrivas/public_data/mvp_20190829/mvp_variant_annots.tsv

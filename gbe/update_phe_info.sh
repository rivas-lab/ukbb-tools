#!/bin/bash

# Finds all `.info` files in the PHE_DIR directory and cats/sorts them together, outputting the catted file to OUT_PATH
# E.g. set PHE_DIR as , and OUT_PATH as `/oak/stanford/groups/mrivas/users/guhan/all_phes.tsv`

PHE_DIR=$1
OUT_PATH=$2

# so that the hash character is interpreted by sort
LANG=C

# finds all the info files, ensures formatting, concatenates them to desired output
find $PHE_DIR -name "*.info" | xargs -i cat {} | awk -F'\t' 'NF == 14' | sort -u > ${OUT_PATH}/phenotype_info.tsv

# also updates icdinfo
awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$8,$2,$8,$8,"Y"}' ${OUT_PATH}/phenotype_info.tsv | tr " " "_" > ${OUT_PATH}/icdinfo.txt

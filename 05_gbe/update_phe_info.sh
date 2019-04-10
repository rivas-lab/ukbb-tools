#!/bin/bash

# Finds all `.info` files in the PHE_DIR directory and cats/sorts them together, outputting the catted file to OUT_PATH
# E.g. set PHE_DIR as , and OUT_PATH as `/oak/stanford/groups/mrivas/users/guhan/all_phes.tsv`

PHE_DIR=$1
OUT_PATH=$2

# so that the hash character is interpreted by sort
LANG=C

# finds all the info files, ensures formatting, concatenates them to desired output
find ${PHE_DIR} -name "*.info" | awk -F'\t' 'NF == 14' | sort -k1,1 -k4,4nr | awk '!_[$1]++'  > ${OUT_PATH}/phenotype_info.tsv

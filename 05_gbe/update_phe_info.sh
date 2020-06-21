#!/bin/bash

# Finds all `.info` files in the PHE_DIR directory and cats/sorts them together, outputting the catted file to OUT_PATH
# E.g. set PHE_DIR as /oak/stanford/groups/mrivas/ukbb24983/phenotypedata, and OUT_DIR as `.`

# sets default directory to $OAK/ukbb24983/phenotypedata, and output to right here
if [ $# -gt 0 ]; then
  PHE_DIR=$1
else
  PHE_DIR="/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/current/info/"
fi

# output redirect
if [ $# -gt 1 ]; then
  OUT_DIR=$2
else 
  OUT_DIR=`pwd`
fi

# so that the hash character is interpreted by sort
LANG=C

# finds all the info files, ensures formatting, concatenates them to desired output
find ${PHE_DIR} -name "*.info" | xargs -i cat {} | awk -F'\t' 'NF == 17' | sort -k1,1 -k4,4nr | awk '!_[$1]++'  > ${OUT_DIR}/phenotype_info.tsv

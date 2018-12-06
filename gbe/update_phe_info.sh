#!/bin/bash

# Finds all `.info` files in the PHE_DIR directory and cats/sorts them together, outputting the catted file to OUT_PATH
# E.g. set PHE_DIR as , and OUT_PATH as `/oak/stanford/groups/mrivas/users/guhan/all_phes.tsv`

PHE_DIR=$1
OUT_PATH=$2

find $PHE_DIR -name "*.info" | xargs -i cat {} | sort > 
#!/bin/bash

# Usage: bash query_master_mrp.sh {mode} {outputprefix} gene1 gene2 ...

# Ex: bash query_master_mrp.sh exome output A1BG WTIP ...

#exome or cal
mode=$1

#output file prefix
out=$2

shift
shift

grepstring=""

for var in "$@"; do
    grepstring+="${var}|"
done

grepstr=${grepstring::-1}

echo $grepstr

cat <(head -1 /oak/stanford/groups/mrivas/ukbb24983/cal/mrp/extras/adjusted_biomarkers/white_british/BIN10030500_gene.tsv) <(egrep -w "$grepstr" /oak/stanford/groups/mrivas/ukbb24983/${mode}/mrp/${mode}_mrp.tsv) >${out}.tsv

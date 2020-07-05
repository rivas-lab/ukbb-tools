#!/bin/bash
set -beEuo pipefail

gwas_current="/oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current"

{

echo "#symlink file"
find ${gwas_current} -name "*INI*logistic*" | sort | while read f ; do echo $f $(readlink -f $f) ; done

} | tr ' ' '\t' > 4_old_files.tsv

find ${gwas_current} -name "*INI*logistic*" | while read f ; do
    rm $(readlink -f $f)
    rm $f
done

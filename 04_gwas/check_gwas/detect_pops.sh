#!/bin/bash

for file in $(ls output | grep check_array_gwas.5); do
    pop=$(tail -n +2 output/$file | head -1 | cut -f3 | cut -d/ -f11)
    filename=$(basename -- "$file")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo $file ${filename}_${pop}.${extension}
done

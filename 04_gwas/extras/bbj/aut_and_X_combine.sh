#!/bin/bash

mkdir -p combined_aut_X
mkdir -p combined_aut_X/bin
mkdir -p combined_aut_X/qt

#biomarkers
#for file in $(ls plink_format/qt/*autosome.txt.gz | grep -v smoking | grep -v cigs | grep -v BMI); do 
#    newfile="${file/plink_format/combined_aut_X}"
#    output="${newfile%.*}"
#    newfile="${file/autosome/combined}"
#    cat <(zcat $file) <(zcat "${file/autosome/chrX}" | tail -n +2) >$newfile
#    gzip $newfile
#done

#other qt
#for file in $(ls plink_format/qt/*autosome*male.txt.gz); do
#    newfile="${file/plink_format/combined_aut_X}"
#    output="${newfile%.*}"
#    newfile="${output/autosome/combined}"
#    cat <(zcat $file) <(zcat "${file/autosome/chrX}" | tail -n +2) >$newfile
#    gzip $newfile
#done

#bins
for file in 

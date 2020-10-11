#!/bin/bash
set -beEuo pipefail

ml load plink2

pfile="/oak/stanford/groups/mrivas/ukbb24983/array-combined/pgen/ukb24983_cal_hla_cnv"
out_d="/oak/stanford/groups/mrivas/ukbb24983/array-combined/afreq_20201009"

pops=('white_british' 'non_british_white' 'african' 's_asian' 'e_asian' 'related' 'others')

################################################
# 500k
################################################

out="${out_d}/$(basename ${pfile})"

plink2 --silent --threads 6 --memory 60000 \
--pfile ${pfile} vzs --out ${out} \
--geno-counts zs cols=chrom,pos,ref,alt,homref,refalt,altxy,hapref,hapalt,missing,nobs

mv ${out}.log ${out}.gcount.log

plink2 --silent --threads 6 --memory 60000 \
--pfile ${pfile} vzs --out ${out} \
--freq zs

mv ${out}.log ${out}.afreq.log

################################################
# Each pop
################################################

for pop in ${pops[@]} ; do

    echo ${pop}

    keep="/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200828/ukb24983_${pop}.phe"
    out="${out_d}/$(basename ${pfile}).${pop}"


    plink2 --silent --threads 6 --memory 60000 \
    --pfile ${pfile} vzs --keep ${keep} --out ${out} \
    --geno-counts zs cols=chrom,pos,ref,alt,homref,refalt,altxy,hapref,hapalt,missing,nobs

    mv ${out}.log ${out}.gcount.log


    plink2 --silent --threads 6 --memory 60000 \
    --pfile ${pfile} vzs --keep ${keep} --out ${out} \
    --freq zs

    mv ${out}.log ${out}.afreq.log

done

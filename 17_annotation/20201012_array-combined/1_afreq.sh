#!/bin/bash
set -beEuo pipefail

ml load plink2

pfile="/oak/stanford/groups/mrivas/ukbb24983/array-combined/pgen/ukb24983_cal_hla_cnv"
sqc_f="/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200828/ukb24983_master_sqc.20200828.phe"
out_d="/oak/stanford/groups/mrivas/ukbb24983/array-combined/annotation/afreq_20201012/plink_output"
arrays=('UKBB' 'UKBL')
pops=('white_british' 'non_british_white' 'african' 's_asian' 'e_asian' 'related' 'others')

plink_wrapper () {
    plink2 --silent --threads 6 --memory 60000 --pfile ${pfile} vzs $@
}

plink_geno_counts () {
    plink_wrapper --geno-counts zs cols=chrom,pos,ref,alt,homref,refalt,altxy,hapref,hapalt,missing,nobs $@
}

plink_freq () {
    plink_wrapper --freq zs $@
}

if [ ! -d ${out_d} ] ; then mkdir -p ${out_d} ; fi

################################################
# 500k
################################################

out="${out_d}/$(basename ${pfile})"

plink_geno_counts --out ${out}
mv ${out}.log ${out}.gcount.log

plink_freq --out ${out}
mv ${out}.log ${out}.afreq.log

################################################
# Each pop
################################################

for pop in ${pops[@]} ; do
    echo ${pop}
    keep="/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200828/ukb24983_${pop}.phe"

    out="${out_d}/$(basename ${pfile}).${pop}"

    plink_geno_counts --keep ${keep} --out ${out}
    mv ${out}.log ${out}.gcount.log

    plink_freq --keep ${keep} --out ${out}
    mv ${out}.log ${out}.afreq.log

done

################################################
# Per UKBL/UKBB array
################################################
for array in ${arrays[@]} ; do
    echo ${array}
    out="${out_d}/$(basename ${pfile}).${array}"

    cat ${sqc_f} | awk -v array=${array} '($8 == array){print $1, $2}' \
    | plink_geno_counts --keep /dev/stdin --out ${out}
    mv ${out}.log ${out}.gcount.log

    cat ${sqc_f} | awk -v array=${array} '($8 == array){print $1, $2}' \
    | plink_freq  --keep /dev/stdin --out ${out}
    mv ${out}.log ${out}.afreq.log

done

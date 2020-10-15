#!/bin/bash
set -beEuo pipefail

ml load plink2

pfile="/oak/stanford/groups/mrivas/ukbb24983/array-combined/pgen/ukb24983_cal_hla_cnv"
out_d="/oak/stanford/groups/mrivas/ukbb24983/array-combined/afreq_20201012/plink_output"
snp_qc_f="/oak/stanford/groups/mrivas/ukbb24983/snp/snp_download/ukb_snp_qc.txt"

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

# one array

for array in 0 1 ; do 
    out="${out_d}/$(basename ${pfile}).array${array}"

    cat ${snp_qc_f} | awk -v array=${array} '($9 == array){print $1}' \
    | plink_geno_counts --out ${out} --extract /dev/stdin

    mv ${out}.log ${out}.gcount.log

    cat ${snp_qc_f} | awk -v array=${array} '($9 == array){print $1}' \
    | plink_freq --out ${out} --extract /dev/stdin

    mv ${out}.log ${out}.afreq.log
done

# both arrays

out="${out_d}/$(basename ${pfile}).both_arrays"

cat ${snp_qc_f} | awk '($9 != 2){print $1}' \
| plink_geno_counts --out ${out} --exclude /dev/stdin

mv ${out}.log ${out}.gcount.log

cat ${snp_qc_f} | awk '($9 != 2){print $1}' \
| plink_freq --out ${out} --exclude /dev/stdin

mv ${out}.log ${out}.afreq.log

################################################
# Each pop
################################################

for pop in ${pops[@]} ; do
    echo ${pop}
    keep="/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200828/ukb24983_${pop}.phe"

    # one array

    for array in 0 1 ; do 
        out="${out_d}/$(basename ${pfile}).${pop}.array${array}"

        cat ${snp_qc_f} | awk -v array=${array} '($9 == array){print $1}' \
        | plink_geno_counts --keep ${keep} --out ${out} --extract /dev/stdin

        mv ${out}.log ${out}.gcount.log

        cat ${snp_qc_f} | awk -v array=${array} '($9 == array){print $1}' \
        | plink_freq --keep ${keep} --out ${out} --extract /dev/stdin

        mv ${out}.log ${out}.afreq.log
    done

    # both arrays

    out="${out_d}/$(basename ${pfile}).${pop}.both_arrays"

    cat ${snp_qc_f} | awk '($9 != 2){print $1}' \
    | plink_geno_counts --keep ${keep} --out ${out} --exclude /dev/stdin

    mv ${out}.log ${out}.gcount.log

    cat ${snp_qc_f} | awk '($9 != 2){print $1}' \
    | plink_freq --keep ${keep} --out ${out} --exclude /dev/stdin

    mv ${out}.log ${out}.afreq.log
done

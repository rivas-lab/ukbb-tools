#!/bin/bash
set -beEuo pipefail

ml load plink2

# input
pfile="/scratch/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/ukb24983_exomeOQFE"
sqc_d="/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200828"

# output
out_d="/scratch/groups/mrivas/ukbb24983/exome/annotation/20201025_exome_oqfe_2020/afreq_hwe"

# parameters
pops=('white_british' 'non_british_white' 'african' 's_asian' 'e_asian' 'related' 'others')

# functions

plink_wrapper () {
    plink2 --silent --threads 6 --memory 60000 --pfile ${pfile} vzs $@
}

plink_geno_counts () {
    plink_wrapper --geno-counts zs cols=chrom,pos,ref,alt,homref,refalt,altxy,hapref,hapalt,missing,nobs $@
}

plink_freq () {
    plink_wrapper --freq zs $@
}

plink_hardy () {
    plink_wrapper --hardy zs midp cols=chrom,pos,ax,gcounts,hetfreq,sexaf,femalep,p $@
}

if [ ! -d ${out_d} ] ; then mkdir -p ${out_d} ; fi

################################################
# 500k
################################################

out="${out_d}/$(basename ${pfile})"

#plink_geno_counts --out ${out}
#mv ${out}.log ${out}.gcount.log

#plink_freq --out ${out}
#mv ${out}.log ${out}.afreq.log

plink_hardy --out ${out}
mv ${out}.log ${out}.hardy.log

################################################
# Each pop
################################################

for pop in ${pops[@]} ; do
    echo ${pop}
    keep="${sqc_d}/ukb24983_${pop}.phe"
    out="${out_d}/$(basename ${pfile}).${pop}"

    plink_geno_counts --keep ${keep} --out ${out}
    mv ${out}.log ${out}.gcount.log

    plink_freq --keep ${keep} --out ${out}
    mv ${out}.log ${out}.afreq.log
        
    plink_hardy --keep ${keep} --out ${out}
    mv ${out}.log ${out}.hardy.log

done


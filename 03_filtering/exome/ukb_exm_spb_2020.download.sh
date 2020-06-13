#!/bin/bash
set -beEuo pipefail 

# download the exome data (2020)
# source data field:
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=23175

log=$(basename $0 .sh).$(date +%Y%m%d-%H%M%S).log

{
    ukbgene evc -m -c1 
    ukbgene evc -c1 
} | tee -a $log

echo $log

# bim file
wget  -nd  biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_spb_exm_chrall_v1.bim

# convert to plink2 pgen/pvar/psam
ml load plink2/20200607
out_prefix=ukb_exm_spb_2020

mv ukb24983_evc_chr1_v1_s49953.fam  ${out_prefix}.fam
mv ukb_evc_chr1_v1.bed              ${out_prefix}.bed
mv ukb_spb_exm_chrall_v1.bim        ${out_prefix}.bim

plink2 --memory 60000 --threads 6 --bfile ${out_prefix} --make-pgen vzs --out ${out_prefix}
mv ${out_prefix}.log ${out_prefix}.pgen.log


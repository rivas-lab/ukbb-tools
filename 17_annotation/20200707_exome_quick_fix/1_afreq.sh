#!/bin/bash
set -beEuo pipefail

ml load plink2/20200620

plink2 --memory 8000 --threads 2 --pfile /oak/stanford/groups/mrivas/ukbb24983/exome/pgen/ukb24983_exome vzs --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200313/ukb24983_white_british.phe --freq cols=chrom,pos,ref,alt,altfreq,nobs zs --out data/ukb24983_exome.white_british

mv data/ukb24983_exome.white_british.log data/ukb24983_exome.white_british.afreq.log

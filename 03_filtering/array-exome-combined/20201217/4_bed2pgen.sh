#!/bin/bash
set -beEuo pipefail

ml load plink2

file=/scratch/groups/mrivas/ukbb24983/array-exome-combined/pgen/ukb24983_cal_hla_cnv_exomeOQFE

plink2 --memory 60000 --threads 8 --bfile ${file} --make-pgen --out ${file} 

mv ${file}.log ${file}.pgen.log

exit 0
## job submission instruction
sbatch -p mrivas --qos=high_p --time=2-0:0:00 --mem=64000 --nodes=1 --cores=8 --job-name=bed2pgen --output=logs/bed2pgen.%A.out --error=logs/bed2pgen.%A.err 4_bed2pgen.sh


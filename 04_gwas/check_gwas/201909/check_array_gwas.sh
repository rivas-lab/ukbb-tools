#!/bin/bash

#SBATCH --job-name=check_array_gwas
#SBATCH   --output=output/check_array_gwas.%A.out
#SBATCH    --error=output/check_array_gwas.%A.err
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH -p owners,normal,mrivas
#SBATCH --nodes=1
#SBATCH --cores=3
#SBATCH --mem=12000
#################
set -beEu -o pipefail

pop=$1

date >&2
bash check_gwas.sh $pop ukb24983_v2_hg19 genotyped
date >&2

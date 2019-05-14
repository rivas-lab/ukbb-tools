#!/bin/bash

#SBATCH --job-name=check_exome_gwas
#SBATCH   --output=check_exome_gwas.%A.out
#SBATCH    --error=check_exome_gwas.%A.err
#SBATCH --time=3:00:00
#SBATCH --qos=normal
#SBATCH -p owners,normal,mrivas
#SBATCH --nodes=1
#SBATCH --cores=3
#SBATCH --mem=12000
#SBATCH --mail-type=END,FAIL
#################
set -beEu -o pipefail

date >&2
bash check_gwas.sh white_british ukb24983_v2_hg38 exome-spb
date >&2


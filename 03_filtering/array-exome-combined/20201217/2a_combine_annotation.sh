#!/bin/bash
set -beEuo pipefail

Rscript 2a_combine_annotation.R

exit 0
#######################
sbatch -p mrivas --nodes=1 --mem=300000 --cores=2 --time=6:00:00 \
    --job-name=combine_annotation --output=logs/combine_annotation.%A.out --error=logs/combine_annotation.%A.err \
    2a_combine_annotation.sh

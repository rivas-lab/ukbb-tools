#!/bin/bash
set -beEuo pipefail

bash $(dirname $0)/check_gwas.sh white_british ukb24983_v2_hg38 exome-spb.PHENO1


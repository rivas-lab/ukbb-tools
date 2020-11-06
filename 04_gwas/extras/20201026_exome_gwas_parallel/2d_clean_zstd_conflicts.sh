#!/bin/bash
set -beEuo pipefail

find /scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE -type f | sort -V > 2d_clean_zstd_conflicts.$(date +%Y%m%d-%H%M%S).txt


#!/bin/bash
set -beEuo pipefail

wc_tsv="gwas-current-gz-wc.$(date +%Y%m%d-%H%M%S).tsv"

{
    echo "#symlink location wc_l" | tr " " "\t"
    find gwas-current-gz-wc -type f -name "*.txt" | sort | parallel cat {}
} > ${wc_tsv}

echo "${wc_tsv}"

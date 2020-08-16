#!/bin/bash
set -beEuo pipefail

ml load R/3.6 gcc

pvar_f="/oak/stanford/groups/mrivas/ukbb24983/array-combined/pgen/ukb24983_cal_hla_cnv.pvar.zst"

master_one_array="/oak/stanford/groups/mrivas/ukbb24983/sqc/one_array_variants.txt"

one_array="$(dirname ${pvar_f})/one_array_variants.txt"
both_arrays="$(dirname ${pvar_f})/both_arrays_variants.txt"

# "one array" file
# generate a symlink
ln -sf ${master_one_array} ${one_array}

# "both arrays" file
Rscript /dev/stdin ${pvar_f} ${one_array} ${both_arrays} << EOF
args <- commandArgs(trailingOnly=TRUE)
suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))

pvar_f <- args[1]
one_array  <- args[2]
both_arrays <- args[3]

one_array_l <- fread(one_array, head=F) %>% pull()

vars <- fread(cmd=paste('zstdcat', pvar_f), select=c('ID')) %>%
mutate(
    is_in_one_array_l = (ID %in% one_array_l)
)

vars %>%
filter(! is_in_one_array_l) %>%
select(ID) %>%
fwrite(both_arrays, col.names=F, row.names=F)
EOF

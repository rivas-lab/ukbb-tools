#!/bin/bash
set -beEuo pipefail

ml load R/3.6 gcc

data_d='/scratch/groups/mrivas/ukbb24983/array-combined/gwas-qc-SE02'
out_f="/oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current/gwas.qc-SE02.$(date +%Y%m%d-%H%M%S).tsv"

Rscript /dev/stdin ${data_d} ${out_f} << EOF
suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))

args <- commandArgs(trailingOnly=TRUE)
data_d <- args[1]
out_f  <- args[2]

pops <- c('white_british', 'non_british_white', 'african', 's_asian', 'e_asian', 'related', 'others', 'metal') 

pops %>% lapply(function(pop){

    lgc_df  <- fread(file.path(data_d, sprintf('%s.qc.tsv', pop)))  %>% rename('GBE_ID'='#GBE_ID')
    cnt_df  <- fread(file.path(data_d, sprintf('%s.cnt.tsv', pop))) %>% rename('GBE_ID'='#GBE_ID')
    icdinfo_df <- fread(sprintf('/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/05_gbe/extras/20200812_GBE_category/icdinfo.20200812/icdinfo.%s.tsv', pop)) %>%
    select(-GBE_short_name_len)

    cnt_df %>%
    left_join(icdinfo_df, by='GBE_ID') %>%
    left_join(
        lgc_df %>%
        mutate(freq_bin = if_else(freq_bin == 'common', 'lgc.common', paste0('lgc', freq_bin))) %>%
        spread(freq_bin, lambda_gc),
        by=c('GBE_ID', 'population')
    )
    
}) %>% bind_rows() -> df

df %>%
rename('#GBE_ID' = 'GBE_ID') %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)
EOF

echo ${out_f}

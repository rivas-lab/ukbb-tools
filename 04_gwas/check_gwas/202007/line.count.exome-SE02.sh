#!/bin/bash
set -beEuo pipefail

ml load R/3.6 gcc

in_f=$1

# example
# in_f='/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/current/white_british/ukb24983_v2_hg38.INI25465.exome-spb.glm.linear.gz'
# out_f='/scratch/groups/mrivas/ukbb24983/exome/gwas-qc/white_british/white_british.INI25465.cnt.tsv'

data_d='/scratch/groups/mrivas/ukbb24983/exome/gwas-qc-SE02'

# configure the output file name
pop=$(basename $(dirname $in_f))
GBE_ID=$(basename $in_f | sed -e 's/ukb24983_v2_hg38.//g' | sed -e 's/.exome-spb//g' | sed -e 's/.glm.logistic.hybrid.gz//g' | sed -e 's/.glm.linear.gz//g' | sed -e 's/.metal.tsv.gz//g')
out_f=${data_d}/${pop}/${pop}.${GBE_ID}.cnt.tsv

# generate the output directory if not exist
if [ ! -d ${data_d}/${pop} ] ; then mkdir -p ${data_d}/${pop} ; fi

# if the output file already exist, exit
if [ -s ${out_f} ] ; then exit 0 ; fi

Rscript /dev/stdin ${out_f} ${in_f} << EOF
suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))
args <- commandArgs(trailingOnly=TRUE)

out_f      <- args[1]
in_f       <- args[2]

# get population and GBE_ID from file path

pop <- basename(dirname(in_f))

GBE_ID <- str_replace_all(
    basename(in_f),
    '^ukb24983_v2_hg38.|.exome-spb|.glm.logistic.hybrid.gz$|.glm.linear.gz$|.metal.tsv.gz$',
    ''
)

# read the files

in_f %>% fread() -> df

# apply SE filter
if('LOG(OR)_SE' %in% colnames(df)){
    df %>% rename('SE'='LOG(OR)_SE') -> df
}
df %>% drop_na(P) %>% filter(SE < .2) -> df

# count lines

data.frame(
    GBE_ID          = GBE_ID,
    population      = pop,
    n_lines         = df %>% nrow(),
    n_non_NA_lines  = df %>% drop_na(P) %>% nrow(),
    n_hits          = df %>% drop_na(P) %>% filter(as.numeric(P) < 5e-8) %>% nrow(),
    stringsAsFactors=F
) -> qc_stats_df

# write to a file
qc_stats_df %>%
rename('#GBE_ID' = 'GBE_ID') %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)

EOF

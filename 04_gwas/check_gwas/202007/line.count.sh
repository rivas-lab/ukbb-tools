#!/bin/bash
set -beEuo pipefail

# ml load R/3.6 gcc

in_f=$1

# example
# in_f='/oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current/white_british/ukb24983_v2_hg19.INI25465.array-combined.glm.linear.gz'
# out_f='/scratch/groups/mrivas/ukbb24983/array-combined/gwas-qc/white_british/white_british.INI25465.cnt.tsv'

data_d='/scratch/groups/mrivas/ukbb24983/array-combined/gwas-qc'

# list of LD independent variants
ld_indep_f="${data_d}/ld_indep.lst"
if [ ! -f ${ld_indep_f} ] ; then
    zcat /oak/stanford/groups/mrivas/private_data/ukbb/variant_filtering/variant_filter_table.20200701.tsv.gz \
    | awk -v FS='\t' '($20=="TRUE"){print $5}' > ${ld_indep_f}
fi

# configure the output file name
pop=$(basename $(dirname $in_f))
GBE_ID=$(basename $in_f | sed -e 's/ukb24983_v2_hg19.//g' | sed -e 's/.array-combined//g' | sed -e 's/.glm.logistic.hybrid.gz//g' | sed -e 's/.glm.linear.gz//g' | sed -e 's/.metal.tsv.gz//g')
out_f=${data_d}/${pop}/${pop}.${GBE_ID}.cnt.tsv

# generate the output directory if not exist
if [ ! -d ${data_d}/${pop} ] ; then mkdir -p ${data_d}/${pop} ; fi

# if the output file already exist, exit
if [ -s ${out_f} ] ; then exit 0 ; fi

Rscript /dev/stdin ${out_f} ${in_f} ${ld_indep_f} << EOF
suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))
args <- commandArgs(trailingOnly=TRUE)

out_f      <- args[1]
in_f       <- args[2]
ld_indep_f <- args[3]

# get population and GBE_ID from file path

pop <- basename(dirname(in_f))

GBE_ID <- str_replace_all(
    basename(in_f),
    '^ukb24983_v2_hg19.|.array-combined|.glm.logistic.hybrid.gz$|.glm.linear.gz$|.metal.tsv.gz$',
    ''
)

# read the files

ld_indep_f %>% fread(header=F) %>% pull() -> ld_indep_lst
in_f %>% fread() -> df

# count lines

data.frame(
    GBE_ID          = GBE_ID,
    population      = pop,
    n_lines         = df %>% nrow(),
    n_non_NA_lines  = df %>% drop_na(P) %>% nrow(),
    n_hits          = df %>% drop_na(P) %>% filter(as.numeric(P) < 5e-8) %>% nrow(),
    n_ld_indep_hits = df %>% drop_na(P) %>% filter(as.numeric(P) < 5e-8, ID %in% ld_indep_lst) %>% nrow(),
    stringsAsFactors=F
) -> qc_stats_df

# write to a file
qc_stats_df %>%
rename('#GBE_ID' = 'GBE_ID') %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)

EOF

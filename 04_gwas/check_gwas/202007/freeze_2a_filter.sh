#!/bin/bash
set -beEuo pipefail

# ml load R/3.6 gcc

in_f=$1

# example
# in_f='/oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current/white_british/ukb24983_v2_hg19.INI25465.array-combined.glm.linear.gz'
# out_f='/scratch/groups/mrivas/ukbb24983/array-combined/gwas/freeze/white_british/white_british.INI25465.1e-3.tsv'

data_d='/scratch/groups/mrivas/ukbb24983/array-combined/gwas/freeze'
p_thr='1e-3'

# configure the output file name
pop=$(basename $(dirname $in_f))
GBE_ID=$(basename $in_f | sed -e 's/ukb24983_v2_hg19.//g' | sed -e 's/.array-combined//g' | sed -e 's/.glm.logistic.hybrid.gz//g' | sed -e 's/.glm.linear.gz//g' | sed -e 's/.metal.tsv.gz//g')
out_f=${data_d}/${pop}/${pop}.${GBE_ID}.${p_thr}.tsv

# generate the output directory if not exist
if [ ! -d ${data_d}/${pop} ] ; then mkdir -p ${data_d}/${pop} ; fi

# if the output file already exist, exit
if [ -s ${out_f} ] ; then exit 0 ; fi

Rscript /dev/stdin ${out_f} ${in_f} ${p_thr} << EOF
suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))
args <- commandArgs(trailingOnly=TRUE)

out_f <- args[1]
in_f  <- args[2]
p_thr <- as.numeric(args[3])

# get population and GBE_ID from file path
pop <- basename(dirname(in_f))

GBE_ID <- str_replace_all(
    basename(in_f),
    '^ukb24983_v2_hg19.|.array-combined|.glm.logistic.hybrid.gz$|.glm.linear.gz$|.metal.tsv.gz$',
    ''
)

# read and filter the summary statistics
df <- fread(
    cmd=paste('zcat', in_f, '|', "sed -e 's/LOG(OR)_SE/SE/g'"),
    select=c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'A1', 'OBS_CT', 'BETA', 'OR', 'SE', 'P', 'ERRCODE'),
    colClasses=c('#CHROM'='character', 'P'='numeric', 'ERRCODE'='character'),
) %>%
rename('CHROM'='#CHROM', 'Variant_ID'='ID') %>%
filter(P <= 1e-3) %>%
mutate(population = pop, GBE_ID = GBE_ID)

if('ERRCODE' %in% colnames(df)){
    df %>% filter(ERRCODE == '.') %>% select(-ERRCODE) -> df
}

if('OR' %in% colnames(df)){
    # convert OR to BETA
    df %>% mutate(BETA = log(as.numeric(OR))) %>% select(-OR) -> df
}

df %>%
select(CHROM, POS, Variant_ID, GBE_ID, population, REF, ALT, A1, OBS_CT, BETA, SE, P) %>%
rename('#CHROM' = 'CHROM') %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)

EOF

echo ${out_f}
exit 0


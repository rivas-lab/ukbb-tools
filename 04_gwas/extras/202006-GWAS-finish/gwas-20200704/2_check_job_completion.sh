#!/bin/bash
set -beEuo pipefail

jobs_f='1_gwas.jobs.tsv'
data_dir='data'
annotated_f="2_gwas.jobs.status.$(date +%Y%m%d-%H%M%S).tsv"
pheno_info_f="/oak/stanford/groups/mrivas/users/$USER/repos/rivas-lab/ukbb-tools/05_gbe/phenotype_info.tsv"

ml load R/3.6 gcc

####################################################################
# R script start
####################################################################
Rscript /dev/stdin ${jobs_f} ${data_dir} ${annotated_f} ${pheno_info_f} << EOF
suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))
####################################################################
args <- commandArgs(trailingOnly=TRUE)
#input
jobs_f       <- args[1]
data_dir     <- args[2]
# output
annotated_f  <- args[3]
pheno_info_f <- args[4]
message(sprintf('jobs_f:       %s', jobs_f))
message(sprintf('data_dir:     %s', data_dir))
message(sprintf('annotated_f:  %s', annotated_f))
message(sprintf('pheno_info_f: %s', pheno_info_f))
####################################################################
# read files
fread(jobs_f) %>%
rename('GBE_ID'='#GBE_ID') -> jobs_df
logs_list     <- fread(cmd=sprintf('find %s -type f -name "*.array-combined.log"',    data_dir), col.names=c('log'))
sumstats_list <- fread(cmd=sprintf('find %s -type f -name "*array-combined.glm*.gz"', data_dir), col.names=c('sumstats'))
pheno_info_f %>% fread(select=c('#GBE_ID', 'PATH')) %>%
rename('GBE_ID'='#GBE_ID') %>%
mutate(
    sumstats_real_dir = str_replace(dirname(dirname(PATH)), '/phenotypedata/', '/array-combined/gwas/')
) %>%
select(-PATH)-> phenotype_df
# merge into a df
jobs_df %>%
mutate(
    log=sprintf('data/%s/ukb24983_v2_hg19.%s.array-combined.log', population, GBE_ID),
    sumstats=sprintf(
        'data/%s/ukb24983_v2_hg19.%s.array-combined.glm.%s.gz', 
        population, GBE_ID,
        if_else(
            str_replace_all(GBE_ID, '[0-9]', '') %in% c('INI', 'QT_FC'),
            'linear', 'logistic.hybrid'
        )
    )
) %>%
left_join(
    logs_list %>% mutate(hasLog=T),
    by='log'
) %>%
left_join(
    sumstats_list %>% mutate(hasSumstats=T),
    by='sumstats'
) %>%
left_join(phenotype_df, by='GBE_ID') %>%
mutate(
    sumstats_real = file.path(sumstats_real_dir, population, basename(sumstats)),
    sumstats_syml = file.path('/oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current', population, basename(sumstats))    
) %>%
select(-sumstats_real_dir) %>%
replace_na(list(hasLog=F, hasSumstats=F)) -> df
# print the counts
df %>% count(hasLog, hasSumstats) %>% print()
# write the results
df %>% 
rename('#GBE_ID' = 'GBE_ID') %>%
fwrite(annotated_f, sep='\t', na = "NA", quote=F)
EOF
####################################################################
# R script end
####################################################################

echo ${annotated_f}

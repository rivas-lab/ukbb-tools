suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

# input
GWAS_covar <- '/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200828/ukb24983_GWAS_covar.20200828.phe'
prev_master_phe <- '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.20200709.phe'
# output
new_master_phe <- '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.20200828.phe'

############################
# main

# read GWAS covars
GWAS_covar %>%
fread(colClasses=c('#FID'='character', 'IID'='character')) %>%
rename('FID'='#FID') -> GWAS_covar_df

# read the older version of master phe file
prev_master_phe %>%
fread(colClasses=c('#FID'='character', 'IID'='character')) %>%
rename('FID'='#FID') -> prev_master_phe_df

# join
GWAS_covar_df %>%
full_join(
    prev_master_phe_df %>%
    select(all_of(c('FID', 'IID', sort(setdiff(colnames(prev_master_phe_df), colnames(GWAS_covar_df)))))),
    by=c('FID', 'IID')
) -> new_master_phe_df

# write to a file
new_master_phe_df %>%
rename('#FID' = 'FID') %>%
fwrite(new_master_phe, sep='\t', na = "NA", quote=F)

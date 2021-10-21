fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

####################################################################
in_f <- args[1]
out_f <- args[2]
split_nonWB_f <- '/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20211020/ukb24983_GWAS_covar.20211020.phe'
####################################################################

split_nonWB_f %>%
fread(select=c('#FID', 'IID', 'split_nonWB'), colClasses = 'character')%>%
rename_with(function(x){str_replace(x, '#', '')}, starts_with("#")) -> split_nonWB_df

in_f %>%
fread(colClasses = c('#FID'='character','FID'='character', 'IID'='character')) %>%
rename_with(function(x){str_replace(x, '#', '')}, starts_with("#")) -> data_df

colname1st <- colnames(data_df)[1]
split_idx <- which(colnames(data_df) == 'split')

if(ncol(data_df) == split_idx){
    out_cols <- c(colnames(data_df), 'split_nonWB')
}else{
    out_cols <- c(colnames(data_df)[1:split_idx], 'split_nonWB', colnames(data_df)[(split_idx+1):ncol(data_df)])
}

data_df %>%
left_join(split_nonWB_df, by=c('FID','IID')) %>%
select(all_of(out_cols)) %>%
rename(!!sprintf('#%s', colname1st) := all_of(colname1st)) %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)

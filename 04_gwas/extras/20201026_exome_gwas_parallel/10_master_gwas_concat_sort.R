fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

####################################################################
out_f <- args[1]
in_f1 <- args[2]

if(length(args)>2){
    in_f2 <- args[3]
}else{
    in_f2 <- NULL # for metal results
}

# out_f <- '/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/ukb24983_exomeOQFE.african.p1e-3.tsv.gz'
# in_f1 <- '/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_filtered_p1e-3/african.linear.tsv.gz'
# in_f2 <- '/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_filtered_p1e-3/african.logistic.hybrid.tsv.gz'

pheno_info_f <- '/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/05_gbe/extras/20200812_GBE_category/GBE_category.20201024.tsv'
var_info_f <- '/oak/stanford/groups/mrivas/ukbb24983/exome/annotation/20201025_exome_oqfe_2020/ukb24983_exomeOQFE.annotation.compact.tsv.gz'

####################################################################

cols_left  <- c('CHROM', 'POS', 'Variant_ID', 'Csq', 'GBE_ID', 'population')
cols_right <- c('SYMBOL', 'GBE_short_name', 'Consequence', 'Gene', 'GBE_category')

####################################################################

message_with_timestamp <- function(msg){
    message(sprintf('[%s] %s', format(Sys.time(), format="%Y%m%d-%H%M%S"), msg))
}

# read files

message_with_timestamp('reading pheno_info_f ..')

pheno_info_f %>%
fread(select=c('#GBE_category', 'GBE_ID', 'GBE_short_name'), colClasses = 'character') %>%
rename('GBE_category'='#GBE_category') -> pheno_info_df

message_with_timestamp('reading var_info_f ..')

var_info_f %>%
fread(select=c('ID', 'Csq', 'Consequence', 'Gene', 'SYMBOL'), colClasses = 'character') %>%
rename('Variant_ID'='ID') -> var_info_df

message_with_timestamp('reading the results file(s) ..')

if(is.null(in_f2)){
    # metal results
    in_f1 %>% fread(colClasses = 'character') %>%
    rename('population'='#population') -> in_df
    
}else{
    in_f1 %>% fread(colClasses = 'character') %>%
    rename('population'='#population') %>%
    rename('T_or_Z_STAT'='T_STAT') -> in_df1
    
    in_f2 %>% fread(colClasses = 'character') %>%
    rename('population'='#population') %>%
    rename('T_or_Z_STAT'='Z_STAT', 'SE'='LOG(OR)_SE') %>%
    mutate(BETA = as.character(log(as.numeric(OR)))) -> in_df2
    
    # each pop. read two files (linear and logistic regression results)
    bind_rows(in_df1, in_df2) -> in_df
    rm(in_df1)
    rm(in_df2)
}

# combine

message_with_timestamp('combining and sorting the data frames ..')

in_df %>%
mutate(in_order = 1:n()) %>%
left_join(var_info_df %>% mutate(var_order=1:n()), by='Variant_ID') %>%
left_join(pheno_info_df, by='GBE_ID') %>%
arrange(var_order, in_order) %>%
select(-in_order, -var_order) -> full_df

message_with_timestamp('writing the results to a file ..')

full_df %>%
select(all_of(c(cols_left, setdiff(colnames(full_df), c(cols_left, cols_right)), cols_right))) %>%
rename('#CHROM' = 'CHROM') %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)

message_with_timestamp('R script finished!')

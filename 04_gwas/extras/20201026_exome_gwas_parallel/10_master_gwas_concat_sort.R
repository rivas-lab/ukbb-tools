fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

####################################################################
out_f <- args[1]
chrom <- args[2]
in_f1 <- args[3]

if(length(args)>3){
    in_f2 <- args[4]
}else{
    in_f2 <- NULL # for metal results
}

nThread <- 6

if (! chrom %in% c(as.character(1:22), 'X', 'Y') ) chrom <- NULL

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

cat_or_zcat <- function(f){
    ifelse(endsWith(f, '.zst'), 'zstdcat', ifelse(endsWith(f, '.gz'), 'zcat', 'cat'))
}

fread_headers <- function(f){
    fread(cmd=paste(
        cat_or_zcat(f), f,
        ' | sed -e "s/^#//g"', ' | head -n1'
    ), colClasses = 'character')
}

fread_w_chr_filter <- function(f, select=NULL, chr=NULL, chrcol=NULL, nThread=1){
    message_with_timestamp(sprintf('  fread_w_chr_filter(chr=%s) ..', chr))
    if(is.null(chrcol)){
        chrcol <- match(c('CHROM'), colnames(fread_headers(f)))
    }
    cmdstr = paste(
        cat_or_zcat(f), f,
        ifelse(
            is.null(chr), '', 
            sprintf(" | awk -v chr=%s -v chrcol=%d '((NR == 1) || ($chrcol == chr))'", chr, chrcol)
        )
    )
    message_with_timestamp(sprintf('  cmd: %s', cmdstr))
    fread(cmd=cmdstr, colClasses = 'character', select=select, nThread=nThread)
}

message_with_timestamp(sprintf('tempdir() = %s', tempdir()))

# read files

message_with_timestamp(sprintf('  in_f1: %s', in_f1))
message_with_timestamp(sprintf('  in_f2: %s', in_f2))

message_with_timestamp('reading pheno_info_f ..')

pheno_info_f %>%
fread(select=c('#GBE_category', 'GBE_ID', 'GBE_short_name'), colClasses = 'character') %>%
rename('GBE_category'='#GBE_category') -> pheno_info_df

message_with_timestamp('reading var_info_f ..')

var_info_f %>%
fread_w_chr_filter(select=c('ID', 'Csq', 'Consequence', 'Gene', 'SYMBOL'), chr=NULL, nThread=nThread) %>%
# fread_w_chr_filter(select=c('ID', 'Csq', 'Consequence', 'Gene', 'SYMBOL'), chr=chrom, nThread=nThread) %>%
rename('Variant_ID'='ID') -> var_info_df

message_with_timestamp(sprintf('  nrow(var_info_df) = %d', nrow(var_info_df)))
message_with_timestamp('reading the results file(s) ..')

if(is.null(in_f2)){
    # metal results
    fread_w_chr_filter(f=in_f1, chr=chrom, chrcol=3, nThread=nThread) %>%
    rename('population'='#population') -> in_df

    message_with_timestamp(sprintf('  nrow(in_df) = %d', nrow(in_df)))    
    
}else{
    # each pop. read two files (linear and logistic regression results)

    fread_w_chr_filter(f=in_f1, chr=chrom, chrcol=3, nThread=nThread) %>%
    rename('population'='#population') %>%
    rename('T_or_Z_STAT'='T_STAT') -> in_df1
    message_with_timestamp(sprintf('  nrow(in_df1) = %d', nrow(in_df1)))

    fread_w_chr_filter(f=in_f2, chr=chrom, chrcol=3, nThread=nThread) %>%
    rename('population'='#population') %>%
    rename('T_or_Z_STAT'='Z_STAT', 'SE'='LOG(OR)_SE') %>%
    mutate(BETA = as.character(log(as.numeric(OR)))) -> in_df2
    message_with_timestamp(sprintf('  nrow(in_df2) = %d', nrow(in_df2)))

    bind_rows(in_df1, in_df2) -> in_df
    rm(in_df1)
    rm(in_df2)
    message_with_timestamp(sprintf('  nrow(in_df)  = %d', nrow(in_df)))
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
fwrite(out_f, sep='\t', na = "NA", quote=F, scipen=100)

message_with_timestamp('R script finished!')

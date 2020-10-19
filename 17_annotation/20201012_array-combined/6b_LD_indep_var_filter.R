fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

####################################################################
csq           <- args[1] # 'pav'
prune_in_csqs <- args[2] # 'ptv' # ptv,pav
out_f         <- args[3] # 'dev.lst'
####################################################################

prune_in <- '/oak/stanford/groups/mrivas/ukbb24983/array-combined/annotation/ld_indep_20201015/plink_output/ukb24983_cal_hla_cnv.white_british.%s.bool.prune.in'
ldmap_f <- '/scratch/groups/mrivas/ukbb24983/cal/ldmap/ldmap_20201018/ukb24983_cal.white_british.ld_map.0.5r2.tsv.gz'
var_QC_f <-  '/oak/stanford/groups/mrivas/ukbb24983/array-combined/annotation/ld_indep_20201015/ukb24983_cal_hla_cnv.var_QCed.tsv.gz'

####################################################################

var_QC_f %>% fread(colClasses = c('#CHROM' = 'character')) %>% 
rename('CHROM'='#CHROM') %>%
select(-FILTER, -geno_data_source) -> var_QC_df

ldmap_f %>% fread(colClasses = c('#CHR_A' = 'character', 'CHR_B' = 'character')) %>%
rename('CHR_A'='#CHR_A') %>%
filter(SNP_A %in% var_QC_df$ID, SNP_B %in% var_QC_df$ID) -> ldmap_df

prune_in_csqs %>% str_split(',') %>% simplify() %>%
lapply(function(csq){
    sprintf(prune_in, csq) %>% fread(head=F)
}) %>% bind_rows() %>% pull() -> prune_in_lst

c(
    ldmap_df %>% filter(SNP_B %in% prune_in_lst) %>% pull(SNP_A),
    ldmap_df %>% filter(SNP_A %in% prune_in_lst) %>% pull(SNP_B)
) %>%
unique() -> remove_lst

var_QC_df %>%
filter(Csq == csq, ! ID %in% remove_lst) %>%
select(ID) %>%
fwrite(out_f, sep='\t', na = "NA", quote=F, col.names = F)

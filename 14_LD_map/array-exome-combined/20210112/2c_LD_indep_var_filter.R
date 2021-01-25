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

if(prune_in_csqs %in% c('none', 'None', 'NONE', 'null', 'Null', 'NULL') ){
    prune_in_csqs <- NULL
}

####################################################################

source('0_parameters.sh')

prune_in <- file.path(
    ldmap_d, 'plink_output', 
    sprintf('%s.white_british.%s.bool.prune.in', basename(pfile), '%s')
)
var_QC_f <- file.path(ldmap_d, 'ukb24983_cal_hla_cnv_exomeOQFE.input.variants.tsv.gz')
ldmap_f <- file.path(ldmap_d, 'ukb24983_cal_hla_exomeOQFE.white_british.ld_map.0.5r2.tsv.gz')

####################################################################

var_QC_f %>% fread(
    select=c('ID', 'geno_data_source', 'ld_indep_array', 'Csq')
) -> var_QC_df

ldmap_f %>% fread(select=c('SNP_A', 'SNP_B', 'R2')) %>%
filter(SNP_A %in% var_QC_df$ID, SNP_B %in% var_QC_df$ID) -> ldmap_df

var_QC_df %>%
drop_na(ld_indep_array) %>%
filter(ld_indep_array) %>%
pull(ID) -> prune_in_array_lst

if(is.null(prune_in_csqs)){
    prune_in_array_lst -> prune_in_lst
}else{
    prune_in_csqs %>% str_split(',') %>% simplify() %>%
    lapply(function(csq){
        sprintf(prune_in, csq) %>% fread(head=F)
    }) %>% bind_rows() %>% pull() -> prune_in_exome_lst
    c(prune_in_array_lst, prune_in_exome_lst) -> prune_in_lst
}

c(
    ldmap_df %>% filter(SNP_B %in% prune_in_lst) %>% pull(SNP_A),
    ldmap_df %>% filter(SNP_A %in% prune_in_lst) %>% pull(SNP_B)
) %>% unique() -> remove_lst

var_QC_df %>% filter(
    geno_data_source == 'exome200k',
    Csq == csq,
    (! ID %in% remove_lst)
) %>% select(ID) %>%
fwrite(out_f, sep='\t', na = "NA", quote=F, col.names = F)

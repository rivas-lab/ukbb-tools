fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

####################################################################

in_f    <- args[1]
patch_f <- args[2]
out_f   <- args[3]

pvar_f  <- args[4]


####################################################################

cat_or_zcat <- function(f){
    ifelse(endsWith(f, '.zst'), 'zstdcat', ifelse(endsWith(f, '.gz'), 'zcat', 'cat'))
}

fread_CHROM <- function(f, select=NULL){
    # fread(cmd=paste(cat_or_zcat(f), f), colClasses = c('#CHROM'='character'), select=select) %>% rename('CHROM'='#CHROM')
    fread(cmd=paste(cat_or_zcat(f), f), colClasses = 'character', select=select) %>% rename('CHROM'='#CHROM')
}

####################################################################

patch_f %>% fread_CHROM() -> patch_df
in_f    %>% fread_CHROM() -> in_df
pvar_f  %>% fread_CHROM() -> pvar_df

pvar_df %>%
left_join(
    bind_rows(
        in_df %>% group_by(ID) %>% slice(1) %>% ungroup(), 
        patch_df
    ) %>% unique(),
    by=c('CHROM', 'POS', 'ID', 'REF', 'ALT')
) %>%
rename('#CHROM' = 'CHROM') %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)

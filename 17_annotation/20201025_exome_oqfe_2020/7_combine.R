suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

# input
pvar_f <- '/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/ukb24983_exomeOQFE.pvar.zst'
data_d <- '/scratch/groups/mrivas/ukbb24983/exome/annotation/20201025_exome_oqfe_2020'
af_hwe_f   <- file.path(data_d, 'ukb24983_exomeOQFE.afreq_hwe.20201025.pvar.zst')
vep_f      <- file.path(data_d, 'UKBexomeOQFE.vep101.tsv.gz')
liftOver_f <- file.path(data_d, 'UKBexomeOQFE.hg19.tsv.gz')
vep_csq_f  <- file.path('..', 'VEP_consequence_group.tsv')

# output
annot_f         <- file.path(data_d, 'ukb24983_exomeOQFE.annotation.20201217.tsv')
annot_compact_f <- file.path(data_d, 'ukb24983_exomeOQFE.annotation.20201217.compact.tsv')

# constants
compact_fields <- c(
    'CHROM', 'POS', 'ID', 'REF', 'ALT', 'FILTER',
    'Allele', 'Csq', 'Consequence', 'SYMBOL', 'Gene', 'HGVSp',
    'f_miss', 'UKB_white_british_hwe_p', 'UKB_white_british_AF', 'UKB_AF',
    'CHROM_hg19', 'POS_hg19', 'REF_hg19', 'ALT_hg19', 'liftOver_unmapped_reason'
)

# functions
cat_or_zcat <- function(f){
    ifelse(endsWith(f, '.zst'), 'zstdcat', ifelse(endsWith(f, '.gz'), 'zcat', 'cat'))
}

fread_CHROM <- function(f, select=NULL){
    fread(
        cmd=paste(cat_or_zcat(f), f, "| sed -e 's/^chr//g'"),
#         cmd=paste(cat_or_zcat(f), f),
        colClasses = c('#CHROM'='character'), select=select
    ) %>% rename('CHROM'='#CHROM')
}

message_with_timestamp <- function(msg){
    message(sprintf('[%s] %s', format(Sys.time(), format="%Y%m%d-%H%M%S"), msg))
}


# read files
message_with_timestamp('Reading the source data files ..')
vep_csq_f %>% fread(select=c('#Consequence', 'Csq')) %>% rename('Consequence'='#Consequence') -> vep_csq_df
pvar_f     %>% fread_CHROM() -> pvar_df
liftOver_f %>% fread_CHROM() -> liftOver_df
af_hwe_f   %>% fread_CHROM() -> af_hwe_df
vep_f      %>% fread_CHROM() -> vep_df
message_with_timestamp('Size of the source dfs are:')
message_with_timestamp(sprintf('  %s: %d x %d', 'pvar_df    ', nrow(pvar_df),     ncol(pvar_df)))
message_with_timestamp(sprintf('  %s: %d x %d', 'vep_df     ', nrow(vep_df),      ncol(vep_df)))
message_with_timestamp(sprintf('  %s: %d x %d', 'af_hwe_df  ', nrow(af_hwe_df),   ncol(af_hwe_df)))
message_with_timestamp(sprintf('  %s: %d x %d', 'liftOver_df', nrow(liftOver_df), ncol(liftOver_df)))

# join
message_with_timestamp('join (1) VEP')
pvar_df %>%
left_join(vep_df, by=c('CHROM', 'POS', 'ID', 'REF', 'ALT')) %>%
left_join(vep_csq_df, by='Consequence') -> full_df
rm(vep_df)
rm(vep_csq_df)
message_with_timestamp(sprintf('  %s: %d x %d', 'full_df    ', nrow(full_df),     ncol(full_df)))

message_with_timestamp('join (2) AF HWE')
full_df %>%
left_join(af_hwe_df, by=c('CHROM', 'POS', 'ID', 'REF', 'ALT'))  -> full_df
rm(af_hwe_df)
message_with_timestamp(sprintf('  %s: %d x %d', 'full_df    ', nrow(full_df),     ncol(full_df)))

message_with_timestamp('join (3) liftOver')
full_df %>%
left_join(liftOver_df, by=c('CHROM', 'POS', 'ID', 'REF', 'ALT')) -> full_df
rm(liftOver_df)
message_with_timestamp(sprintf('  %s: %d x %d', 'full_df    ', nrow(full_df),     ncol(full_df)))

# write the results
message_with_timestamp('writing the results to disk ..')
full_df %>%
select(all_of(c(compact_fields, setdiff(colnames(full_df), compact_fields)))) %>%
rename('#CHROM' = 'CHROM') %>%
fwrite(annot_f, sep='\t', na = "NA", quote=F)
message(annot_f)

full_df %>%
select(all_of(compact_fields)) %>%
rename('#CHROM' = 'CHROM') %>%
fwrite(annot_compact_f, sep='\t', na = "NA", quote=F)
message(annot_compact_f)

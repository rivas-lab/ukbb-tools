suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

# filenames
ukb_d <- '/scratch/groups/mrivas/ukbb24983'
# output
var_annot_f <- file.path(ukb_d, 'array-exome-combined/annotation/20210111/ukb24983_cal_hla_cnv_exomeOQFE.annot_20210111.tsv')
var_annot_compact_f <- file.path(dirname(var_annot_f), 'ukb24983_cal_hla_cnv_exomeOQFE.annot_compact_20210111.tsv')
# input
array_f <- file.path(ukb_d, 'array-combined/annotation/annotation_20201012/ukb24983_cal_hla_cnv.annot_20201023.tsv.gz')
exome_f <- file.path(ukb_d, 'exome/annotation/20201025_exome_oqfe_2020/ukb24983_exomeOQFE.annotation.20210108.tsv.gz')
pvar_f  <- file.path(ukb_d, 'array-exome-combined/pgen/ukb24983_cal_hla_cnv_exomeOQFE.pvar.gz')
ldindep <- file.path(ukb_d, 'array-exome-combined/ldmap/ldmap_20210112/ukb24983_cal_hla_exomeOQFE.white_british.bool.in.gz')

# constants
cols <- c(
    'CHROM', 'POS', 'ID', 'REF', 'ALT', 
    'FILTER', 'POS_total', 'Allele', 'Csq', 'Consequence', 
    'SYMBOL', 'Gene', 'ld_indep', 'geno_data_source', 'array',
    'CNV_POS_s', 'CNV_POS_e', 'UKB_white_british_MAF', 
    'UKB_white_british_hwe_p', 'mgi_notes', 'f_miss', 'f_miss_UKBB', 'f_miss_UKBL',
    'LoF', 'LoF_filter', 'LoF_flags', 'LoF_info'
)

# function
get_POS_total <- function(df){
    # compute the location on the linear coordinate system (chr1-22, X, Y, MT) for plotting
    # original: https://github.com/rivas-lab/ukbb-tools/blob/master/17_annotation/20201012_array-combined/7_finalize.ipynb
    df %>% select(CHROM, POS, ID) %>%
    mutate(CHROM_X = if_else(CHROM == 'XY', 'X', CHROM)) -> CHROM_POS_df

    CHROM_POS_df %>% group_by(CHROM_X) %>%
    summarise(chr_len = max(POS), .groups = 'drop') %>%
    left_join(data.frame(CHROM_X = c(1:22, 'X', 'Y', 'MT'), CHROM_order=1:25, stringsAsFactors=F), by='CHROM_X') %>%
    arrange(CHROM_order) %>% mutate(CHROM_tot= cumsum(as.numeric(chr_len)) - chr_len) %>%
    select(CHROM_X, CHROM_tot) %>% left_join(CHROM_POS_df, by='CHROM_X') %>%
    mutate(POS_total = POS + CHROM_tot) %>%
    select(ID, POS_total)    
}

show_log <- function(msg){
    message(sprintf('[%s] %s', format(Sys.time(), format="%Y%m%d-%H%M%S"), msg))
}

# read LD indep results
ldindep %>% fread(head=F) %>% pull() -> ld_indep_vars
show_log(sprintf('%d LD indep variants', length(ld_indep_vars)))

# read pvar file
pvar_f %>%
fread(colClasses = c('#CHROM'='character', 'ID'='character')) %>%
rename('CHROM'='#CHROM') %>%
mutate(sort_order = 1:n()) -> pvar_df

pvar_df %>% filter(geno_data_source != 'exome200k') %>%
pull(ID) -> ID_array

pvar_df %>% filter(geno_data_source == 'exome200k') %>%
pull(ID) -> ID_exome

show_log(sprintf('%d array + %d exome variants are detected in the pvar file', length(ID_array), length(ID_exome)))

# read & filter array
array_f %>% fread(colClasses = 'character', sep='\t') %>%
rename('CHROM'='#CHROM') %>% filter(ID %in% ID_array)  %>%
rename('UKB_white_british_hwe_p' = 'hwe_p') %>% mutate(
    UKB_white_british_MAF = as.numeric(UKB_white_british_MAF)
) %>% select(-geno_data_source, -POS_total, -ld_indep, -FILTER) -> array_df

# read & filter exome
exome_f %>% fread(colClasses = 'character', sep='\t') %>%
rename('CHROM'='#CHROM') %>% filter(ID %in% ID_exome) %>% rename(
    'CHROM_hg38'='CHROM', 'POS_hg38'='POS', 'REF_hg38'='REF', 'ALT_hg38'='ALT',
    'CHROM'='CHROM_hg19', 'POS'='POS_hg19', 'REF'='REF_hg19', 'ALT'='ALT_hg19'
) %>% mutate(
    UKB_white_british_MAF = pmin(1 - as.numeric(UKB_white_british_AF), as.numeric(UKB_white_british_AF))
) %>% select(-liftOver_unmapped_reason, -FILTER) -> exome_df

# combine
bind_rows(array_df, exome_df) %>% select(-CHROM, -POS, -REF, -ALT) %>%
inner_join(pvar_df, by='ID') %>%
mutate(ld_indep = (ID %in% ld_indep_vars)) -> combined_df

show_log(sprintf('combined_df of %d x %d', nrow(combined_df), ncol(combined_df)))

# sort
combined_df %>% left_join(get_POS_total(combined_df), by='ID') %>%
arrange(sort_order) %>% select(-sort_order) -> full_df

show_log(sprintf('full_df of %d x %d', nrow(full_df), ncol(full_df)))

# write
full_df %>% select(all_of(c(cols, 'HGVSp'))) %>%
rename('#CHROM' = 'CHROM') %>%
fwrite(var_annot_compact_f, sep='\t', na = "NA", quote=F)

full_df %>% select(all_of(c(cols, setdiff(colnames(full_df), cols)))) %>%
rename('#CHROM' = 'CHROM') %>%
fwrite(var_annot_f, sep='\t', na = "NA", quote=F)


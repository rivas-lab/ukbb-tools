fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

####################################################################
# annot_f <- '/oak/stanford/groups/mrivas/ukbb24983/exome/annotation/20201025_exome_oqfe_2020/ukb24983_exomeOQFE.annotation.20201217.compact.tsv.gz'
in_annot_f  <- args[1]
out_annot_f <- args[2]
dup_f   <- '/oak/stanford/groups/mrivas/ukbb24983/exome/qc/oqfe_2020/intermediate_files/ukb24983_exomeOQFE.duplicates.tsv.gz'

####################################################################

dup_f %>% fread(colClasses = c('#CHROM'='character')) %>%
rename('CHROM'='#CHROM') -> dup_df

in_annot_f %>%
fread(colClasses = c('#CHROM'='character')) %>%
rename('CHROM'='#CHROM') %>%
mutate(
    QC_duplicated = ! (ID %in% (dup_df$ID)),
    QC_sample_miss = (f_miss < .1),
    QC_WB_HWE_p = (log10(UKB_white_british_hwe_p)>-15)
) %>%
mutate(FILTER = paste0(
    FILTER,
    if_else(QC_duplicated  & (!is.na(QC_duplicated)),  '', ';duplicated'),
    if_else(QC_sample_miss & (!is.na(QC_sample_miss)), '', ';sample_miss'),
    if_else(QC_WB_HWE_p    & (!is.na(QC_WB_HWE_p)),    '', ';WB_HWE_p')
)) %>%
mutate(FILTER = str_replace(FILTER, '^.;', '')) %>%
select(
    - QC_duplicated, - QC_sample_miss, - QC_WB_HWE_p
) -> annot_df

annot_df %>%
rename('#CHROM' = 'CHROM') %>%
fwrite(out_annot_f, sep='\t', na = "NA", quote=F)
message(out_annot_f)

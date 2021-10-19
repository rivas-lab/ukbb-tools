# this script is adapted based on 
# ../array-combined/array-combined.p_factor.v5.R

suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))


# input
Csq_f <- '/oak/stanford/groups/mrivas/ukbb24983/cal/annotation_20200912/ukb24983_cal_cALL_v2_hg19.vep101.noLoF.Csq.tsv.gz'
HLA_f <- '/oak/stanford/groups/mrivas/ukbb24983/hla/pgen/ukb_hla_v3.pvar'
CNV_f <- '/oak/stanford/groups/mrivas/ukbb24983/cnv/pgen/cnv.pvar'
imp_f_template <- '/oak/stanford/groups/mrivas/ukbb24983/imp/pgen/var_QC_noHWE/ukb24983_imp_chr%s_v3_QC_noHWE.pvar.zst'
clinvar_f <- '../array-combined/clinvar_20200914_patho.tsv'
# this is generated with clinvar_extract.sh

# output
p_factor <- '/oak/stanford/groups/mrivas/ukbb24983/array_imp_combined/snpnet/penalty.v1.rds'

##################
# main

clinvar_df <- fread(clinvar_f) %>%
rename_with(function(x){str_replace(x, '#', '')}, starts_with("#")) %>%
mutate(CLNSIG_curated = str_replace_all(CLNSIG, '[,/].+', '')) %>%
select(-CLNSIG)

HLA_df <- fread(HLA_f, select=c('ID', 'ALT')) %>%
mutate(ID_ALT = paste(ID, ALT, sep='_'), w=.75)

CNV_df <- fread(CNV_f, select=c('ID', 'ALT')) %>%
mutate(ID_ALT = paste(ID, ALT, sep='_'), w=1)

imp_df <- c(1:22, 'X', 'XY') %>%
lapply(function(chr){
    imp_f <- sprintf(imp_f_template, chr)
    fread(
        cmd = sprintf("zstdcat %s | grep -vE '^##'", imp_f), 
        select=c('ID', 'ALT')
    )
}) %>% bind_rows() %>%
mutate(ID_ALT = paste(ID, ALT, sep='_'), w=1)

Csq <- fread(Csq_f, select=c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'Csq')) %>%
rename_with(function(x){str_replace(x, '#', '')}, starts_with("#"))

Csq %>%
mutate(CHROM_X = if_else(CHROM == 'XY', 'X', CHROM)) %>%
left_join(clinvar_df, by=c('CHROM_X'='CHROM', 'POS', 'REF', 'ALT')) %>%
select(-CHROM_X) %>%
replace_na(list(CLNSIG_curated='')) %>%
mutate(
    w = if_else(
        Csq == 'ptv' | CLNSIG_curated == 'Pathogenic', .5,
        if_else(
            Csq == 'pav' | CLNSIG_curated == 'Likely_pathogenic', .75, 1
        )
    )
) -> Csq_w

bind_rows(
    Csq_w %>%
    mutate(ID_ALT = paste(ID, ALT, sep='_')) %>%
    select(ID_ALT, w),
    
    HLA_df %>% select(ID_ALT, w),
    
    CNV_df %>% select(ID_ALT, w),
    
    imp_df %>% select(ID_ALT, w)
) -> weights

weights %>%
deframe() %>%
saveRDS(file = p_factor)

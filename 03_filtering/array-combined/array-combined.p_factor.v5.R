suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))


# input
Csq_f <- '/oak/stanford/groups/mrivas/ukbb24983/cal/annotation_20200912/ukb24983_cal_cALL_v2_hg19.vep101.noLoF.Csq.tsv.gz'
HLA_f <- '/oak/stanford/groups/mrivas/ukbb24983/hla/pgen/ukb_hla_v3.pvar'
CNV_f <- '/oak/stanford/groups/mrivas/ukbb24983/cnv/pgen/cnv.pvar'
clinvar_f <- 'clinvar_20200914_patho.tsv'
# this is generated with clinvar_extract.sh

# output
p_factor <- '/oak/stanford/groups/mrivas/ukbb24983/array-combined/snpnet/penalty.v5.rds'

##################
# main

HLA_df <- fread(HLA_f) %>%
rename('CHROM'='#CHROM') %>%
mutate(ID_ALT = paste(ID, ALT, sep='_'), w=.75)

CNV_df <- fread(CNV_f) %>%
rename('CHROM'='#CHROM') %>%
mutate(ID_ALT = paste(ID, ALT, sep='_'), w=1)

Csq <- fread(Csq_f) %>%
rename('CHROM'='#CHROM')

clinvar_df <- fread(clinvar_f) %>% rename('CHROM'='#CHROM') %>%
mutate(CLNSIG_curated = str_replace_all(CLNSIG, '[,/].+', ''))

Csq %>%
mutate(CHROM_X = if_else(CHROM == 'XY', 'X', CHROM)) %>%
left_join(clinvar_df%>% select(-CLNSIG), by=c('CHROM_X'='CHROM', 'POS', 'REF', 'ALT')) %>%
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
    
    CNV_df %>% select(ID_ALT, w)
) -> weights

weights %>%
deframe() %>%
saveRDS(file = p_factor)


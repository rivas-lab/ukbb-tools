suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

# output
p_factor_v2 <- '/oak/stanford/groups/mrivas/ukbb24983/array-combined/snpnet/penalty.v2.rds'
# input
p_factor_v1 <- '/oak/stanford/groups/mrivas/ukbb24983/array-combined/snpnet/penalty.rds'
hla_pvar    <- '/oak/stanford/groups/mrivas/ukbb24983/hla/pgen/ukb_hla_v3.p.pvar'

# get the list of HLA allelotypes from pvar file
hla_pvar %>% fread() %>% rename('CHROM'='#CHROM') %>%
mutate(ID_ALT = paste(ID, ALT, sep='_')) %>% pull(ID_ALT) -> hla_ID_ALT

# read the p.factor rds file (v1)
p_factor_v1 %>% readRDS() %>% enframe() -> p_factor_v1_df

print('p.factor (v1)')
p_factor_v1_df %>% count(value) %>% print()

print('p.factor (v1) HLA allelotype')
p_factor_v1_df %>% filter(name %in% hla_ID_ALT) %>% count(value) %>% print()

# assign .75 (coding weights) to HLA allelotypes
p_factor_v1_df %>% 
mutate(value = if_else(name %in% hla_ID_ALT, .75, value)) -> p_factor_v2_df

print('p.factor (v2)')
p_factor_v2_df %>% count(value) %>% print()

# save the new data
p_factor_v2_df %>% deframe() %>% saveRDS(file = p_factor_v2)

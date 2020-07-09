suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

tab_f <-  '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/2005693/37855/download/ukb37855.tab'
cols_f <- paste0(tab_f, '.columns')

out_f <- 'data/ukb.2005693.37855.coding339.long.tsv'

GBE_IDs <- c(
    'INI21048','INI21049','INI21051','INI21052',
    'INI21053','INI21054','INI21055','INI21056',
    'INI21058','INI21059','INI21060','INI21061'
)

GBE_IDs %>% 
lapply(function(x){as.integer(str_replace(x, 'INI', ''))}) -> UKB_Fields

cols_df <- fread(cols_f, skip=1, head=F)

cols_df %>% filter(V3 %in% UKB_Fields) -> cols_sub_df

tab_df <- fread(
    tab_f, 
    select=c('f.eid', cols_sub_df$V1),
    colClasses=c('f.eid'='character')
)

tab_df %>%
gather(field, val, -f.eid) %>%
drop_na(val) %>%
rename('#FID' = 'f.eid') %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)

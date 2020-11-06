suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))
args <- commandArgs(trailingOnly=TRUE)

in_f <- args[1]
out_f <- args[2]

in_f %>%
fread(colClasses = c('#file'='character')) %>%
rename('file'='#file') %>%
separate(file, c(rep(NA, 9), 'population', 'basename'), remove=F, sep='/') %>%
separate(basename, c(NA, 'GBE_ID'), remove=F, fill='left', sep='\\.', extra='drop') %>%
mutate(GBE_CAT = str_replace_all(GBE_ID, '[0-9]', '')) %>%
select(population, GBE_ID, GBE_CAT, NF_wc_l, wc_l, non_NA_wc_l, basename, file) %>%
rename('#population' = 'population') %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)

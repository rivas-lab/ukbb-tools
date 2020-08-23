suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))

args <- commandArgs(trailingOnly=TRUE)

# input
long_df_f <- args[1]
coding_f  <- args[2]
# output
long_df_annot_f <- args[3]

# read HC coding file 
## cf. https://github.com/rivas-lab/ukbb-phenotyping/blob/master/README.md

coding_f %>% fread(sep='|', col.names=c('GBE_NAME', 'coding6', 'coding19', 'HC_number')) %>%
mutate(GBE_ID = paste0('HC', HC_number), coding6 = as.character(coding6)) %>%
select(-HC_number) %>%
mutate(
    coding19 = str_replace_all(coding19, '\\.', '')
) -> coding

bind_rows(
    coding %>%
    select(GBE_ID, coding6) %>%
    separate(coding6, paste0('coding6_', 1:400), sep=',', remove=T, fill='right', extra='drop') %>%
    gather(tmp, val, -GBE_ID) %>%
    filter(!is.na(val)) %>%
    select(-tmp) %>%
    mutate(coding = 6),
    
    coding %>%
    select(GBE_ID, coding19) %>%
    separate(coding19, paste0('coding19_', 1:400), sep=',', remove=T, fill='right', extra='drop') %>%
    gather(tmp, val, -GBE_ID) %>%
    filter(!is.na(val)) %>%
    select(-tmp) %>%
    mutate(coding = 19)
) %>%
arrange(GBE_ID, coding, val) %>%
select(GBE_ID, coding, val) -> coding_long

# read the extracted file
fread(cmd=paste('zstdcat', long_df_f), colClasses=c('#IID'='character')) %>%
rename('IID'='#IID') -> long_df

# add GBE_ID field, sort, and write to a file
long_df %>%
mutate(coding = if_else(field == 20002, 6, 19)) %>%
left_join(coding_long, by = c("val", "coding")) %>%
arrange(GBE_ID, field, coding, val, IID, time, array, tab_file) %>%
select( GBE_ID, field, coding, val, IID, time, array, tab_file) %>%
rename('#GBE_ID' = 'GBE_ID') %>%
fwrite(long_df_annot_f, sep='\t', na = "NA", quote=F)

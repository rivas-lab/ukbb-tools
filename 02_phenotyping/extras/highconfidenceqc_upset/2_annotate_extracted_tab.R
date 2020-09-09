suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))

args <- commandArgs(trailingOnly=TRUE)

# input
long_df_f <- args[1]
coding_f  <- args[2]
# output
long_df_annot_f <- args[3]

############################
# functions
############################

find_full_ICD_codes <- function(GBEID, coding, ICD10s){
    # In the mapping file, we have set of ICD codes. Some of them refers to the code itself and its descendants.
    # This function generates a regex pattern (ptn) from the information in the mapping file and
    # query the pattern for the all ICD codes (`str_detect`) and report all the ICD10 codes that matches the pattern.
    paste0('^', coding %>% filter(GBE_ID == GBEID) %>% pull(coding19) %>% str_replace_all(',', '|^')) -> ptn
    ICD10s[str_detect(ICD10s, ptn)]
}

self_reported_field_IDs <- c(20002)

############################
# main
############################

# read HC coding file
## cf. https://github.com/rivas-lab/ukbb-phenotyping/blob/master/README.md

coding_f %>% fread(sep='|', col.names=c('GBE_NAME', 'coding6', 'coding19', 'HC_number')) %>%
mutate(GBE_ID = paste0('HC', HC_number), coding6 = as.character(coding6)) %>%
select(-HC_number) %>%
mutate(
    coding19 = str_replace_all(coding19, '\\.', '')
) -> coding

# read the extracted file
fread(cmd=paste('zstdcat', long_df_f), colClasses=c('#IID'='character')) %>%
rename('IID'='#IID') -> long_df

# get the set of ICD10 codes used in our observation
long_df %>% filter(! field %in% self_reported_field_IDs) %>%
pull(val) %>% unique() %>% sort() -> ICD10s

bind_rows(
    # self-reported codes
    coding %>%
    select(GBE_ID, coding6) %>%
    separate(coding6, paste0('coding6_', 1:10), sep=',', remove=T, fill='right', extra='drop') %>%
    gather(tmp, val, -GBE_ID) %>%
    filter(!is.na(val)) %>%
    select(-tmp) %>%
    mutate(coding = 6),

    # ICD codes
    coding %>% pull(GBE_ID) %>% unique() %>% sort() %>%
    lapply(function(GBEID){
        vals = find_full_ICD_codes(GBEID, coding, ICD10s)
        data.frame(
            GBE_ID = rep(GBEID, length(vals)),
            coding = rep(19, length(vals)),
            val = vals,
            stringsAsFactors=F
        )
    }) %>%
    bind_rows()
) %>%
arrange(GBE_ID, coding, val) %>%
select(GBE_ID, coding, val) -> coding_long

# add GBE_ID field, sort, and write to a file
long_df %>%
mutate(coding = if_else(field %in% self_reported_field_IDs, 6, 19)) %>%
left_join(coding_long, by = c("val", "coding")) %>%
arrange(GBE_ID, field, coding, val, IID, time, array, tab_file) %>%
select( GBE_ID, field, coding, val, IID, time, array, tab_file) %>%
rename('#GBE_ID' = 'GBE_ID') %>%
fwrite(long_df_annot_f, sep='\t', na = "NA", quote=F)

fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

####################################################################
# source(file.path(dirname(script.name), 'misc.R'))
####################################################################

out_f <- '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/additional_medications/misc/ukb2007183_ukb40831.field20003.tsv'

tab_and_cols <- setNames(
    list(        
        c(20003)
    ),
    c(
        '/scratch/groups/mrivas/ukbb24983/phenotypedata/2007183/40831/download/ukb40831.tab'
    )
)

####################################################################

read_tab_as_long <- function(tab_file, field_IDs){
    tabc <- fread(
        file=sprintf('%s.columns', tab_file), 
        skip=1, header=F
    )
    colnames(tabc)<-c('colname', 'f', 'field', 'time', 'array')
    
    select_cols <- tabc %>% 
    filter(field %in% field_IDs) %>%
    select(colname) %>% 
    pull()

    fread(
        tab_file, select=c('f.eid', select_cols),
        colClasses=c('f.eid'='character')
    ) %>%
    rename('IID'='f.eid') %>%
    gather(key, val, -IID) %>%
    drop_na(val) %>%
    mutate(key=str_replace(key, 'f.', '')) %>%
    separate(key, c('field', 'time', 'array'))
}

long_df <- tab_and_cols %>%
names() %>%
lapply(function(tab_file){
    cols <- tab_and_cols[[tab_file]]
    read_tab_as_long(tab_file, cols) %>%
    mutate(tab_file=basename(tab_file))
}) %>%
bind_rows()

long_df %>% 
rename('#IID' = 'IID') %>%
fwrite(str_replace_all(out_f, '.zst$', ''), sep='\t', na = "NA", quote=F)

system(paste(
    'zstd', '--rm', '-19',
    str_replace_all(out_f, '.zst$', ''),
    sep=' '
), intern=F, wait=T)

print(sprintf('The extracted phenotype is saved to: %s.zst', str_replace_all(out_f, '.zst$', '')))

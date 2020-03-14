fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

# data_dir <- '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/self_reported_ethnicity'
# in_f     <- file.path(data_dir, 'misc', 'ukb9796_ukb24611_f21000.tsv')
# out_f    <- file.path(data_dir, 'phe',  'ukb9796_ukb24611_f21000.phe')

in_f  <- args[1]
out_f <- args[2]

####################################################################


find_consistent_answers <- function(df){
    # For self-reported answers, find the consistent answers
    df %>% 
    gather("col", "val", -IID) %>%
    drop_na(val) %>%
    group_by(IID) %>%
    summarise(
        vals = list(unique(val)),
        n_vals = length(unique(val)),
        fst = first(unique(val))
    ) %>%
    ungroup() %>%
    filter(n_vals == 1) %>%
    mutate(phe = fst) %>% 
    select(IID, phe)
}

in_f %>% 
fread(header=T, data.table=F) %>%
find_consistent_answers() %>%
rename(f21000 = phe) %>% 
mutate(FID=IID) %>%
select(FID, IID, f21000) %>%
rename('#FID' = 'FID') %>%
fwrite(out_f, sep='\t')

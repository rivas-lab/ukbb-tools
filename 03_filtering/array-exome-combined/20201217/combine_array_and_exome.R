fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

####################################################################
out_f   <- args[1]
array_f <- args[2]
exome_f <- args[3]
pvar_f  <- args[4]
####################################################################

# read pvar file
pvar_f %>%
fread(colClasses = c('#CHROM'='character')) %>%
rename('CHROM'='#CHROM') %>%
mutate(sort_order = 1:n()) -> pvar_df

pvar_df %>% filter(geno_data_source != 'exome200k') %>%
pull(ID) -> ID_array

pvar_df %>% filter(geno_data_source == 'exome200k') %>%
pull(ID) -> ID_exome

# read & filter array
array_f %>%
fread(colClasses = 'character') %>%
rename('CHROM'='#CHROM') %>%
filter(ID %in% ID_array) -> array_df

# read & filter exome
exome_f %>%
fread(colClasses = 'character') %>%
rename('CHROM'='#CHROM') %>%
filter(ID %in% ID_exome) -> exome_df

# combine
pvar_df %>% select(CHROM, POS, ID, REF, ALT, sort_order) %>%
inner_join(
    bind_rows(
        array_df %>% select(-CHROM, -POS, -REF, -ALT), 
        exome_df %>% select(-CHROM, -POS, -REF, -ALT)
    ),
    by='ID'
) %>%
arrange(sort_order) %>%
select(-sort_order) -> full_df

# write the results
full_df %>%
rename('#CHROM' = 'CHROM') %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)

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
fread(
    colClasses = c('#CHROM'='character'),
    select=c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'FILTER', 'geno_data_source')
) %>%
rename('CHROM'='#CHROM') %>%
filter(FILTER == '.') %>%
mutate(sort_order = 1:n()) -> pvar_df

pvar_df %>% filter(geno_data_source != 'exome200k') %>%
pull(ID) -> ID_array

pvar_df %>% filter(geno_data_source == 'exome200k') %>%
pull(ID) -> ID_exome

# read & filter array
array_f %>%
fread(colClasses = 'character') %>%
rename('CHROM'='#CHROM') %>%
filter(ID %in% ID_array) %>%
select(-CHROM, -POS, -REF, -ALT) -> array_df

# read & filter exome
exome_f %>%
fread(colClasses = 'character') %>%
rename('CHROM'='#CHROM') %>%
filter(ID %in% ID_exome) %>%
select(-CHROM, -POS, -REF, -ALT) -> exome_df

# combine
rm(ID_array)
rm(ID_exome)
bind_rows(array_df, exome_df) -> combined_df
rm(array_df)
rm(exome_df)

pvar_df %>% select(CHROM, POS, ID, REF, ALT, sort_order) %>%
inner_join(combined_df, by='ID') %>%
arrange(sort_order) %>% select(-sort_order) -> full_df
rm(combined_df)

# write the results
full_df %>%
rename('#CHROM' = 'CHROM') %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)


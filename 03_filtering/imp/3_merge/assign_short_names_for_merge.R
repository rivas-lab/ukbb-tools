args = commandArgs(trailingOnly=TRUE)

require(tidyverse)
require(data.table)

in_file  <- args[1] # pvar file
prefix   <- args[2]
out_file <- args[3] # tsv file

fread(in_file) %>% 
rename('CHROM' = '#CHROM', full_ID = ID) %>% 
mutate(ID = paste0(prefix, 1:n())) %>%
select(CHROM, POS, full_ID, REF, ALT, ID) %>%
fwrite(out_file, sep='\t')

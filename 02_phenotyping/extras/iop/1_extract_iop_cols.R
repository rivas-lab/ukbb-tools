suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(data.table))

phe_path <- '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata'
source_tab_file <- file.path(phe_path, '2005693/37855/download/ukb37855.tab')
out_tab_file <- file.path(phe_path, 'extras/iop/misc/ukb2005693_ukb37855_IOP.tsv')
GBE_IDs <- c(5254, 5255, 5262, 5263)

selectCols <- c(simplify2array(c(lapply(GBE_IDs, function(x){c(paste0('f.', x, '.0.0'), paste0('f.', x, '.1.0'))}))))
colClases <- setNames(rep('double', length(selectCols)), selectCols)
colClases[['f.eid']] <- 'character'

tab_df <- fread(
    source_tab_file,
    colClasses = colClases,
    select = c('f.eid', selectCols)
) %>%
rename('IID' = 'f.eid')

tab_df %>% 
fwrite(out_tab_file, sep='\t', na = "NA", quote = F)

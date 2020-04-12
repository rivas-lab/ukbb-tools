require(tidyverse)
require(data.table)
require(readxl)

names_df <- as.data.frame(read_excel('GBE_names.xlsx'))
print(file)
pops <- c(
'white_british',
'non_british_white',
'african',
's_asian',
'e_asian'
)

for(pop in pops){
    fread(paste0('icdinfo.', pop, '.txt'), head=F) %>%
    rename(GBE_ID = V1) %>%
    left_join(
        names_df %>% select(-GBE_N, -GBE_category, -GBE_NAME), by='GBE_ID'
    ) %>% 
    fwrite(paste0('icdinfo.', pop, '.shortnames.txt'), sep='\t', col.names=F)
}

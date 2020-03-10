# Script that reshapes OPCS4 codes (f.41272*) and matching event (f.41282*) to a table with columns FID, A012, A013, A014..... [all OPCS]  
suppressMessages(require("tidyverse"))
suppressMessages(require("data.table"))
suppressMessages(require("dplyr"))
suppressMessages(require("tidyr"))
suppressMessages(require("lubridate"))
suppressMessages(require("reshape"))
suppressMessages(require("reshape2"))

### Read in all OPCS4 codes and matcing dates of event
opcs_f <- '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/2005693/37855/download/ukb37855.tab'
opcs_df <- fread(
    cmd = paste0('cat ', opcs_f, ' | cut -f1,12814-12930,13207-13323'), sep='\t',data.table=T
    )
# Select columns with OPCS4
df <- opcs_df %>% select(f.eid, f.41272.0.0:f.41272.0.116)

df_opcs <- df %>% gather(column, OPCS, f.41272.0.0:f.41272.0.116) %>%
    mutate(id = row_number())

df1 <- opcs_df %>% select(f.eid, f.41282.0.0:f.41282.0.116)

# Select columns with dates
df_date <- df1 %>% gather(column, date, f.41282.0.0:f.41282.0.116) %>%
    mutate(id = row_number())

# Join 
full_df <- left_join(df_opcs, df_date, by = c("f.eid" = "f.eid", "id" = "id"))

# Rename and select columns FID, OPCS, date
opcs_df <- full_df %>% 
dplyr::rename('FID' = 'f.eid') %>% 
select(FID, OPCS, date) 

l.sort <- opcs_df[order(opcs_df$FID),]

# Reshape 
print("Reshape dataframe...")
reshape_df <- reshape(l.sort, timevar = "OPCS", idvar = c("FID"), direction = "wide")

# remove "age_event." from all columns
for ( col in 1:ncol(reshape_df)){
    colnames(reshape_df)[col] <-  sub("date.", "", colnames(reshape_df)[col])
}

# Sort columns 
print("Sorted columns")
sorted_df <- reshape_df %>% select(FID, order(names(.)))

# Write table 
sorted_df %>% fwrite(
    '/oak/stanford/groups/mrivas/users/justesen/projects/disease_progression/phe/opcs.formatting.phe',
    sep='\t', na = "NA", quote=F
)
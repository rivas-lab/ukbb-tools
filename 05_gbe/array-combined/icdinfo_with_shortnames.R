require(tidyverse)
require(data.table)
require(readxl)
require(writexl)

names_df <- as.data.frame(read_excel('GBE_names.xlsx', sheet="GBE_names"))

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

df <- fread(cmd=paste0('cat ', 'icdinfo.txt', ' | cut -f1-3'), sep='\t', data.table=F)
colnames(df) <- c('GBE_ID', 'GBE_N', 'GBE_NAME')

df_clean <- df %>% mutate(
    GBE_short_name = GBE_NAME,
    GBE_short_name = str_replace_all(GBE_short_name, '_', ' '),
    GBE_short_name = str_replace_all(GBE_short_name, '\\(left\\)', '(L)'),
    GBE_short_name = str_replace_all(GBE_short_name, '\\(right\\)', '(R)'),
    GBE_short_name = str_replace_all(GBE_short_name, 'percentage', '%'),
    GBE_short_name = str_replace_all(GBE_short_name, '90th percentile', '90%ile'),
    GBE_short_name = str_replace_all(GBE_short_name, 'number', '#'),
    GBE_short_name = str_replace_all(GBE_short_name, 'Number', '#'),
    GBE_short_name = str_replace_all(GBE_short_name, 'Frequency', 'Freq.'),        
    GBE_short_name = str_replace_all(GBE_short_name, 'higher than', '>'),
    GBE_short_name = str_replace_all(GBE_short_name, 'lower than', '<'),    
    GBE_short_name = str_replace_all(GBE_short_name, 'Volume of', 'Vol. of'),
    GBE_short_name = str_replace_all(GBE_short_name, 'volume of', 'vol. of'),    
    GBE_short_name = str_replace_all(GBE_short_name, 'predicted', 'pred.'),
    GBE_short_name = str_replace_all(GBE_short_name, 'blood pressure', 'BP'),
    GBE_short_name = str_replace_all(GBE_short_name, 'Average', 'Ave.'),
    GBE_short_name = str_replace_all(GBE_short_name, 'average', 'ave.'),
    GBE_short_name = str_replace_all(GBE_short_name, 'distance', 'dist.'),
    GBE_short_name = str_replace_all(GBE_short_name, ', automated reading', ' (AR)'),   
    GBE_short_name = str_replace_all(GBE_short_name, 'cholelithiasis/gall stones', 'Gallstones'),
    GBE_short_name = str_replace_all(GBE_short_name, 'Body mass index \\(BMI\\)', 'BMI'),
    GBE_short_name = str_replace_all(GBE_short_name, 'Weighted-mean', 'WA'), # WA. weighted average
    GBE_short_name = str_replace_all(GBE_short_name, 'treatments/medications', 'medications'),
    GBE_short_name = str_replace_all(GBE_short_name, 'Peak expiratory flow \\(PEF\\)', 'PEF'),
    GBE_short_name = str_replace_all(GBE_short_name, 'Forced expiratory volume in 1-second \\(FEV1\\)', 'FEV1'),
    GBE_short_name = str_replace_all(GBE_short_name, 'Forced vital capacity \\(FVC\\)', 'FVC'),
    GBE_short_name = str_replace_all(GBE_short_name, 'statistic', 'stat.'),
    GBE_short_name = str_replace_all(GBE_short_name, "Alzheimer's disease/dementia", "Alzheimer's/dementia"),
    GBE_short_name = str_replace_all(GBE_short_name, 'Time spent outdoors in', 'Outdoor time,'),
    GBE_short_name = str_replace_all(GBE_short_name, 'Nitrogen dioxide', 'NO2'),
    GBE_short_name = str_replace_all(GBE_short_name, 'Particulate matter air pollution', 'air pollution'),
    GBE_short_name = str_replace_all(GBE_short_name, 'sound level of noise pollution', 'noise lvl.'),
    GBE_short_name = str_replace_all(GBE_short_name, 'platelet \\(thrombocyte\\)', 'platelet'),
    GBE_short_name = str_replace_all(GBE_short_name, 'White blood cell \\(leukocyte\\)', 'White blood cell'),
    GBE_short_name = str_replace_all(GBE_short_name, 'Red blood cell \\(erythrocyte\\)', 'Red blood cell'),
    GBE_short_name = str_replace_all(GBE_short_name, 'Age at menopause \\(last menstrual period\\)', 'Age at menopause'),
    GBE_short_name = str_replace_all(GBE_short_name, 'difficulty/problems', 'problems'),
    GBE_short_name = str_replace_all(GBE_short_name, 'Nucleated red blood cell', 'Nuc. red blood cell'),
    GBE_short_name = str_replace_all(GBE_short_name, 'night-time', 'nighttime'),
    GBE_short_name = str_replace_all(GBE_short_name, 'air pollution', 'air poll.'),
    GBE_short_name = str_replace_all(GBE_short_name, ';', ''),       
    GBE_short_name = str_replace_all(GBE_short_name, 'heart attack/myocardial infarction', 'heart attack (MI)'),
    GBE_short_name = str_replace_all(GBE_short_name, 'Childhood sunburn occasions', 'Childhood sunburn'),
    GBE_short_name = str_replace_all(GBE_short_name, 'deep venous thrombosis \\(dvt\\)', 'DVT'),
    GBE_short_name = str_replace_all(GBE_short_name, 'dvt', 'DVT'),
    GBE_short_name = str_replace_all(GBE_short_name, 'Attention deficit \\(hyperactivity\\) disorder \\(ADD/ADHD\\)', 'ADD/ADHD'),
    GBE_short_name = str_replace_all(GBE_short_name, 'COPD \\(chronic obstructive pulmonary disease\\)', 'COPD'),    
    GBE_short_name = str_replace_all(GBE_short_name, 'pulmonary embolism', 'PE'),
    GBE_short_name = str_replace_all(GBE_short_name, 'hereditary/genetic', 'genetic'),
    GBE_short_name = str_replace_all(GBE_short_name, 'C reactive protein', 'C-reactive protein'),
    GBE_short_name = str_replace_all(GBE_short_name, '\\s+$', ''),
    GBE_short_name = str_replace_all(GBE_short_name, '^\\s+', ''),   
    GBE_short_name = firstup(GBE_short_name) # capitalize the first letter
) %>% select(GBE_ID, GBE_N, GBE_NAME, GBE_short_name)

joined <- df_clean %>% full_join(names_df, by=c('GBE_ID','GBE_NAME')) %>%
select(GBE_ID, GBE_N.x, GBE_N.y, GBE_NAME, GBE_short_name.x, GBE_short_name.y, Units_of_measurement) %>% 
mutate(
    GBE_short_name = ifelse(is.na(GBE_short_name.y), GBE_short_name.x, GBE_short_name.y),
    GBE_N = ifelse(is.na(GBE_N.x), GBE_N.y, GBE_N.x),
    GBE_short_name_len = str_length(GBE_short_name),
    GBE_category = str_replace_all(GBE_ID, '[0-9]+', '')
) %>% 
select(GBE_category, GBE_ID, GBE_N, GBE_NAME, GBE_short_name, GBE_short_name_len, Units_of_measurement)

new_df_clean <- joined %>% filter(GBE_ID %in% df_clean$GBE_ID)

new_df_clean %>% fwrite(
    'icdinfo.shortnames.tsv', sep='\t'
)

write_xlsx(joined, "GBE_names.xlsx")

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

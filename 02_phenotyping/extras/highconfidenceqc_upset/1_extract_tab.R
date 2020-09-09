fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

####################################################################
# source(file.path(dirname(script.name), 'misc.R'))
####################################################################

out_f <- '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/highconfidenceqc_upset/ukb37855_ukb40831_icd.tsv'

# ukb40831
#
# - 20002: Non-cancer illness code, self-reported
# - 41202: Diagnoses - main ICD10
# - 41204: Diagnoses - secondary ICD10
# - 40001: Underlying (primary) cause of death: ICD10
# - 40002: Contributory (secondary) causes of death: ICD10
# - 41201: External causes - ICD10
#
# ukb37855
# - 41270: Diagnoses - ICD10

tab_and_cols <- setNames(
    list(
        c(20002, 41202, 41204, 40001, 40002, 41201),
        c(41270)
    ),
    c(
        '/scratch/groups/mrivas/ukbb24983/phenotypedata/2007183/40831/download/ukb40831.tab',
        '/scratch/groups/mrivas/ukbb24983/phenotypedata/2005693/37855/download/ukb37855.tab'
#         '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/2007183/40831/download/ukb40831.tab',
#         '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/2005693/37855/download/ukb37855.tab'
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
fwrite(out_f, sep='\t', na = "NA", quote=F)

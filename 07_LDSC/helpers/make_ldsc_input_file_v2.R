fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

####################################################################
# source(file.path(dirname(script.name), 'misc.R'))
####################################################################
read_and_filter_plink_sumstats <- function(in_f){
    df1 <- fread(in_f, colClasses=c('#CHROM'='character', 'ID'='character', 'P'='character')) %>%
    rename('CHROM'='#CHROM')%>%
    mutate(A2 = if_else(A1 == ALT, REF, ALT)) %>%
    filter(TEST=='ADD', REF != 'N', ALT != 'N')%>%
    drop_na('P')

    rename_dict <- list()
    rename_dict[['LOG(OR)_SE']] <- 'SE'
    df_new_colnames <- df1 %>% colnames() %>%
    lapply(function(x){ifelse(x %in% names(rename_dict), rename_dict[[x]], x)})
    colnames(df1) <- df_new_colnames

    if(!'BETA' %in% colnames(df1)){
        df2 <- df1 %>% drop_na(OR)
        df2$BETA <- log(as.numeric(df2$OR))
    }else{
        df2 <- df1 %>% drop_na(BETA)
    }
    df2 %>% drop_na(BETA)
}

read_ldscore_files <- function(ldscore_d){
    Sys.glob(file.path(ldscore_d, "*.ldscore.gz")) %>% 
    lapply(function(ldscore_f){
        ldscore_df <- fread(ldscore_f, colClasses=c('CHR'='character'))
        ldscore_df %>%
        left_join(
            ldscore_df %>% 
            count(CHR, BP) %>% 
            filter(n == 1) %>%
            select(CHR, BP),
            by=c('CHR', 'BP')
        )
    }) %>%
    bind_rows()
}

####################################################################
out_f      <- args[1]
in_f       <- args[2]
ldscore_d  <- args[3]

####################################################################

ldscore_df <- read_ldscore_files(ldscore_d)

df <- read_and_filter_plink_sumstats(in_f)

merged_df <- df %>%
select(-ID)%>%
inner_join(
    ldscore_df %>% select(CHR, BP, SNP),
    by=c('CHROM'='CHR', 'POS'='BP')
)%>%
rename('ID'='SNP') %>%
# mutate(P = if_else(P < .Machine$double.xmin, as.character(.Machine$double.xmin), P)) %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)

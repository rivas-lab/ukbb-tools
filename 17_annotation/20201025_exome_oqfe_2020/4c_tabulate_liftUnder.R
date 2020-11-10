suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

# input
data_d <- '/scratch/groups/mrivas/ukbb24983/exome/annotation/20201025_exome_oqfe_2020'
pvar_f <- file.path(data_d, 'UKBexomeOQFE.pvar.zst')
liftOver_mapped_f   <- file.path(data_d, 'UKBexomeOQFE.hg19.mapped.tsv.gz')
liftOver_unmapped_f <- file.path(data_d, 'UKBexomeOQFE.hg19.unmapped.txt.gz')

# output
liftOver_f <- file.path(data_d, 'UKBexomeOQFE.hg19.tsv')

# functions
cat_or_zcat <- function(f){
    ifelse(endsWith(f, '.zst'), 'zstdcat', ifelse(endsWith(f, '.gz'), 'zcat', 'cat'))
}

fread_CHROM <- function(f, select=NULL){
    fread(cmd=paste(cat_or_zcat(f), f), colClasses = c('#CHROM'='character'), select=select) %>% rename('CHROM'='#CHROM')
}

# read unmapped
bind_cols(
    fread(cmd=paste(
        'zcat', liftOver_unmapped_f, '|', "egrep -v '#'", '|', 
        'sed -e "s/^#//g"', '|', 'cut -f4'
    ), head=F, sep='\t') %>% rename('ID'=1),
    
    fread(cmd=paste(
        'zcat', liftOver_unmapped_f, '|', "egrep '#'", '|', 
        'sed -e "s/^#//g"'
    ), head=F, sep='\t') %>% rename('liftOver_unmapped_reason'=1)
) -> liftOver_unmapped_df

# read pvar
pvar_f %>% fread_CHROM() -> pvar_df

# read mapped 
liftOver_mapped_f %>% fread_CHROM() %>%
rename('CHROM_hg19'='CHROM', 'POS_hg19'='POS', 'REF_hg19'='REF', 'ALT_hg19'='ALT') -> liftOver_mapped_df


# join
pvar_df %>%
left_join(bind_rows(
    liftOver_mapped_df,
    liftOver_unmapped_df
), by='ID') -> full_df

# save it to a file
full_df %>%
rename('#CHROM' = 'CHROM') %>%
fwrite(liftOver_f, sep='\t', na = "NA", quote=F)

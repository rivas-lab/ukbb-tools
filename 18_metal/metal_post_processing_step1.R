fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

####################################################################
in_files  <- args[1]
metal_f   <- args[2]
out_f     <- args[3]
####################################################################

# read metal output
metal <- fread(metal_f, colClasses='character')

in_files %>% fread(head=F) %>% pull() %>%
lapply(function(f){
    f %>%
    fread(
        select=c('#CHROM', 'POS', 'ID', 'OBS_CT'),
        colClasses=c('#CHROM'='character', 'POS'='numeric', 'ID'='character', 'OBS_CT'='numeric')
    ) %>%
    rename('CHROM'='#CHROM')    
}) %>%
bind_rows() %>%
group_by(CHROM, POS, ID) %>%
summarise(OBS_CT=sum(OBS_CT)) %>%
ungroup() %>%
as.data.frame() %>%
inner_join(metal %>% rename('ID' = 'MarkerName'), by='ID') %>%
mutate(
    Allele1 = toupper(Allele1),
    Allele2 = toupper(Allele2)
) %>%
rename(
    'A1'  = 'Allele1', 
    'REF' = 'Allele2'
    # Allele1 - the first allele for this marker in the first file where it occurs
    # Allele2 - the second allele for this marker in the first file where it occurs
    #
    # .. However, it seems like the ordering is pretty much random 
    # (we need to apply flipfix)
) %>%
mutate(ALT = A1) %>%
arrange(CHROM, as.numeric(POS)) %>%
select(c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'A1', 'OBS_CT', colnames(metal)[4:ncol(metal)])) %>%
rename('#CHROM' = 'CHROM', 'BETA' = 'Effect', 'SE' = 'StdErr', 'P' = 'P-value') %>%
fwrite(out_f, sep='\t')

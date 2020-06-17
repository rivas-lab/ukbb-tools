fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

####################################################################
in_f  <- args[1]
out_f <- args[2]
finngen_master_bim_f <- '/scratch/groups/mrivas/public_data/summary_stats/finngen_r3/finngen_r3_variants.master.tsv.gz'
####################################################################

in_f %>%
fread(colClasses=c('#chrom'='character')) %>%
rename('chrom'='#chrom') %>%
mutate(hg38_ID=paste(chrom, pos, ref, alt, sep=':')) %>%
select(-chrom, -pos, -ref, -alt) -> in_df

finngen_master_bim_f %>%
fread(
    select=c('#CHROM', 'POS', 'REF', 'ALT', 'ID', 'hg19_REF_ALT_flip', 'hg38_ID'),
    colClasses=c('#CHROM'='character', 'ID'='character', 'hg38_ID'='character', 'hg19_REF_ALT_flip'='numeric')
) %>%
rename('CHROM'='#CHROM') %>%
mutate(ID = if_else(is.na(ID), paste('hg38', hg38_ID, sep=':'), ID)) %>%
right_join(in_df, by='hg38_ID') %>%
filter(hg19_REF_ALT_flip != 0) %>%
mutate( # apply flipfix
    maf          = if_else(hg19_REF_ALT_flip == 1, maf,          1 - maf),
    maf_controls = if_else(hg19_REF_ALT_flip == 1, maf_controls, 1 - maf_controls),
    maf_cases    = if_else(hg19_REF_ALT_flip == 1, maf_cases,    1 - maf_cases),
    beta         = hg19_REF_ALT_flip * beta
) %>%
select(-hg19_REF_ALT_flip) %>%
mutate(OR = exp(beta)) %>%
rename('P'='pval', 'SE'='sebeta') %>%
select(CHROM, POS, REF, ALT, ID, P, OR, SE, hg38_ID, maf, maf_cases, maf_controls, beta, rsids, nearest_genes) %>%
rename('#CHROM' = 'CHROM') %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)

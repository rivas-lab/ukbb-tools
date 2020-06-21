suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

# input

finngen_dir <- '/scratch/groups/mrivas/public_data/summary_stats/finngen_r3'
finngen_vars_f      <- file.path(finngen_dir, 'finngen_r3_variants.tsv.gz')
finngen_vars_hg19_f <- file.path(finngen_dir, 'finngen_r3_variants.hg19.fasta.tsv.gz')
ukb_vars_f <- '/oak/stanford/groups/mrivas/ukbb24983/array_imp_combined/pgen_v2/ukb24983_cal_hla_cnv_imp.pvar.zst'

# output
out_f               <- file.path(finngen_dir, 'finngen_r3_variants.master.tsv')

#################

ukb_vars_df <- fread(cmd=paste('zstdcat', ukb_vars_f), colClasses=c('#CHROM'='character')) %>%
rename('CHROM'='#CHROM')

finngen_vars_hg19_df <- fread(finngen_vars_hg19_f, select=c('#CHROM', 'POS', 'ID', 'FASTA_REF'), colClasses=c('#CHROM'='character', 'ID'='character', 'FASTA_REF'='character')) %>%
rename('CHROM'='#CHROM')

finngen_vars_f %>%
fread(colClasses=c('#CHROM'='character', 'ID'='character')) %>%
rename('hg38_CHROM'='#CHROM', 'hg38_POS'='POS') %>%
left_join(finngen_vars_hg19_df, by='ID') %>%
mutate(hg19_REF_ALT_flip = if_else(REF == toupper(FASTA_REF), 1, if_else(ALT == toupper(FASTA_REF), -1, 0))) %>%
replace_na(list(hg19_REF_ALT_flip=0)) %>%
select(-FASTA_REF) %>%
rename('hg38_REF'='REF', 'hg38_ALT'='ALT', 'hg38_ID'='ID') %>%
mutate(
    REF=if_else(hg19_REF_ALT_flip == 1, hg38_REF, if_else(hg19_REF_ALT_flip == -1, hg38_ALT, 'NA')),
    ALT=if_else(hg19_REF_ALT_flip == 1, hg38_ALT, if_else(hg19_REF_ALT_flip == -1, hg38_REF, 'NA'))
) %>%
left_join(ukb_vars_df, by=c('CHROM','POS','REF','ALT')) %>%
select(CHROM, POS, REF, ALT, ID, hg19_REF_ALT_flip, hg38_CHROM, hg38_POS, hg38_REF, hg38_ALT, hg38_ID, rsids, nearest_genes) %>%
arrange(CHROM, POS, ID) %>%
rename('#CHROM' = 'CHROM') %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)

message(out_f)

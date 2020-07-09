suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

# input
afreq = 'data/ukb24983_exome.white_british.afreq.zst'
pvar = '/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/ukb24983_exome.pvar.zst'
annotation = '/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/spb/data/ukb_exm_spb-white_british_variant_annots.tsv.gz'

# output
out = '/oak/stanford/groups/mrivas/ukbb24983/exome/annotation/ukb_exm_spb_2020.white_british.tsv'

pvar_df <- fread(cmd=paste('zstdcat', pvar))
message('pvar')
annotation_df <- fread(annotation, select=c('CHROM', 'POS', 'REF', 'ALT', 'ID', 'Gene', 'Gene_symbol', 'Consequence'))
message('annoation')
afreq_df <- fread(cmd=paste('zstdcat', afreq), select=c('ID', 'ALT_FREQS'))
message('afreq')

pvar_df %>%
rename('CHROM'='#CHROM') %>%
mutate(pvar_order=1:n()) %>%
left_join(
    afreq_df, by=c('ID')
) %>%
left_join(
    annotation_df %>% rename('rsID'='ID'),
    by=c('CHROM', 'POS', 'REF', 'ALT')
) %>%
arrange(pvar_order) %>%
select(-pvar_order) %>%
rename('#CHROM' = 'CHROM', 'ALT_FREQ_white_british'='ALT_FREQS') %>%
fwrite(out, sep='\t', na = "NA", quote=F)

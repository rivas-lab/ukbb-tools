############################################
# GWAS Manhattan plot
#
#   Yosuke Tanigawa
#   2020/9/23
############################################

fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

source(file.path(dirname(script.name), 'gwas_plot_misc.R'))
####################################################################
gwas_sumstats_f  <- args[1]
annotation_f     <- args[2]
manhattan_plot_f <- args[3]
####################################################################

annot <- read_annotation_tbl(annotation_f)

gwas_sumstats_f %>%
fread(colClasses = c('#CHROM'='character')) %>% rename('CHROM'='#CHROM') %>%
select(CHROM, POS, ID, REF, ALT, P) %>%
mutate(CHROM = if_else(CHROM == 'XY', 'X', CHROM), P = as.numeric(P), POS = as.numeric(POS), log10P = -log10(P)) %>%
left_join(annot %>% select(ID, Gene_symbol), by='ID') %>%
mutate(
    rankP = rank(P),
    repel_label = if_else((rankP <= 30) & (P <= 1e-6), Gene_symbol, ''),
    color = if_else(CHROM %in% c(2 * 1:11, 'Y'), "2_even_chrs", "1_odd_chrs"),
) -> gwas_compact_df

gwas_compact_df %>% compute_gwas_plot_df() %>%
filter(P < 1e-2) %>% plot_manhattan() -> p_manhattan

ggsave(manhattan_plot_f, p_manhattan)

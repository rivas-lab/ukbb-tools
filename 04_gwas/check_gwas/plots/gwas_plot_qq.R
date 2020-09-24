############################################
# GWAS qq plot
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
gwas_sumstats_f <- args[1]
qq_plot_f       <- args[2]
####################################################################

gwas_sumstats_f %>%
fread(colClasses = c('#CHROM'='character')) %>% rename('CHROM'='#CHROM') %>%
qq_plot() -> p_qq

ggsave(qq_plot_f, p_qq)

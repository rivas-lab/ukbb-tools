suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))
# Allele strings are too long in CNV dataset. Let's clean them up.

args <- commandArgs(trailingOnly=TRUE)
in_f <- args[1]
out_f <- args[2]

fread(cmd=paste('cat', in_f, '|', 'cut -f1,2,5-')) %>%
mutate(Allele = if_else(str_length(Allele)>1, '+', '-')) %>%
mutate(HGVSc = str_replace_all(HGVSc, 'ins[ACGT]+', 'ins')) %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)

fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

####################################################################
in_f  <- args[1]
out_f <- args[2]
if(length(args)>2){
    p_thr <- as.numeric(args[3])
}else{
    p_thr <- 1e-3
}
####################################################################
# this script apply P-value threshold (on log10(P) scale to handle P < 1e-300) and save the results to the specified file

in_f %>%
fread(colClasses = c('#CHROM'='character', 'P'='character')) %>%
rename('CHROM'='#CHROM') %>%
separate(P, c('P_base', 'P_exp'), sep='e', remove=F, fill='right') %>%
replace_na(list(P_exp='0')) %>%
mutate(log10P = log10(as.numeric(P_base)) + as.numeric(P_exp)) %>%
filter(-log10P > -log10(p_thr)) %>%
select(-P_base, -P_exp) %>%
rename('#CHROM' = 'CHROM') %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)

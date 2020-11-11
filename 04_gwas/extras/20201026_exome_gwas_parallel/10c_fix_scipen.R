fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

options(scipen=999)

####################################################################
in_f  <- args[1]
out_f <- args[2]
####################################################################

in_f %>%
fread(colClasses = 'character', nThread=2) %>%
mutate(POS = format(as.integer(as.numeric(POS)), scientific = F)) %>%
fwrite(out_f, sep='\t', na = "NA", quote=F, scipen=999, nThread=2)

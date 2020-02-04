fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

####################################################################
# source(file.path(dirname(script.name), 'misc.R'))
####################################################################

in_f  <- args[1]
out_f <- args[2]

fread(cmd=paste('zcat', in_f)) %>% 
mutate(
    FASTA_REF = toupper(FASTA_REF),
    BETA = as.numeric(BETA),
    is_flip = ! (FASTA_REF == REF),
    ALT = if_else(is_flip, REF, ALT),
    REF = if_else(is_flip, A1,  REF),
    A1  = if_else(is_flip, ALT, A1),
    BETA = if_else(is_flip, -1 * BETA, BETA),
    Direction = if_else(
        is_flip, 
        str_replace_all(str_replace_all(str_replace_all(Direction, '-', 'm'), '\\+', '-'), 'm', '+'),
        Direction
    )
) %>%
select(-FASTA_REF, -is_flip) %>%
fwrite(out_f, sep='\t')

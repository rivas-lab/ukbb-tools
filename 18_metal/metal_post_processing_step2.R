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

HLA_allelotypes_with_P_as_ref <- c('DPA1_103', 'DRB3_9901', 'DRB4_9901', 'DRB5_9901')
# For all CNV alleles and most HLA allelotypes, REF == 'N'
# However, the 4 alleloeyptes specified here, REF == 'P'

fread(cmd=paste('zcat', in_f)) %>% 
mutate(
    FASTA_REF = toupper(FASTA_REF),
    BETA = as.numeric(BETA),
    is_flip = if_else(
        # special treatment for 4 HLA allelotypes
        (ID %in% HLA_allelotypes_with_P_as_ref), REF != 'P',
        if_else(
            # special treatment for CNV and HLA alleles
            REF != 'N' & ALT != 'N',
            FASTA_REF != REF,
            REF != 'N'
        )
    ),
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

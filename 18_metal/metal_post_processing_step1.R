fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

####################################################################
# source(file.path(dirname(script.name), 'misc.R'))
####################################################################

loci_f  <- args[1]
metal_f <- args[2]
out_f   <- args[3]

# read loci file
loci <- fread(cmd=paste('zcat', loci_f)) %>%
rename('CHROM' = '#CHROM')

# read metal output
metal <- fread(metal_f)

# join them together and change col names
df <- loci %>% 
mutate(order = 1:n()) %>%
inner_join(
    metal %>% rename('ID' = 'MarkerName'),
    by='ID'
) %>%
mutate(
    Allele1 = toupper(Allele1),
    Allele2 = toupper(Allele2)
) %>%
rename(
    'A1'  = 'Allele1', 
    'REF' = 'Allele2'
    # Allele1 - the first allele for this marker in the first file where it occurs
    # Allele2 - the second allele for this marker in the first file where it occurs
    #
    # .. However, it seems like the ordering is pretty much random 
    # (we need to apply flipfix)
) %>%
mutate(ALT = A1) %>%
arrange(order) %>%
select(c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'A1', colnames(metal)[4:ncol(metal)])) %>%
rename('#CHROM' = 'CHROM', 'BETA' = 'Effect', 'SE' = 'StdErr', 'P' = 'P-value')

# write the results to file
df %>% fwrite(out_f, sep='\t')

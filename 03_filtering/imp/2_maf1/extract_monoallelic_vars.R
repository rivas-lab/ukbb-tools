args = commandArgs(trailingOnly=TRUE)

require(tidyverse)
require(data.table)

in_file  <- args[1]
out_file <- args[2]

#in_file<-'/oak/stanford/groups/mrivas/ukbb24983/imp/pgen/maf1/debug/debug_XY.pvar.zst'
#out_file<-'/oak/stanford/groups/mrivas/ukbb24983/imp/pgen/maf1/debug/debug_XY.lst'


df <- fread(cmd=paste('zstdcat', in_file, sep=' ')) 
uniq_POSs <- df %>% count(POS) %>% filter(n == 1) %>% select(POS) %>% pull()
df %>%filter(POS %in% uniq_POSs) %>% select(ID) %>% fwrite(out_file, col.names=F)

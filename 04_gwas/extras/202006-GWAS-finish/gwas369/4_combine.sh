#!/bin/bash
set -beEuo pipefail
ml load R/3.6 gcc

f1=$1 # the results will be overwritten to f1
f2=$2 # additional GWAS results

############################################################
# tmp dir
############################################################
tmp_dir_root="$LOCAL_SCRATCH"
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
# echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT

############################################################
# filenames
############################################################

tmpf1=${tmp_dir}/$(basename $f1)

cp $f1 $tmpf1

############################################################
# combine the sumstats files
############################################################
Rscript /dev/stdin $tmpf1 $f2 ${f1%.gz} << EOF
args <- commandArgs(trailingOnly=TRUE)
suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))
in_f1 <- args[1]
in_f2 <- args[2]
out_f <- args[3]

df <- bind_rows(
    fread(in_f1, colClasses='character'),
    fread(in_f2, colClasses='character')
) %>%
rename('CHROM'='#CHROM') %>%
mutate(CHROM =  factor(CHROM, levels = c(1:22, 'X', 'XY', 'Y', 'MT'))) %>%
arrange(CHROM, as.integer(POS), ID) %>%
mutate(CHROM =  as.character(CHROM)) %>%
rename('#CHROM' = 'CHROM') %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)
EOF

bgzip -l9 -f ${f1%.gz}

echo ${f1}

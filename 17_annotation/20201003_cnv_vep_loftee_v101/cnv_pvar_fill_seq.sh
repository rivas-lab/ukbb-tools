#!/bin/bash
set -beEuo pipefail

# ref_fa='/oak/stanford/groups/mrivas/public_data/vep/20200912/homo_sapiens/101_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa'
ref_fa='/scratch/groups/mrivas/public_data/vep/20200912/homo_sapiens/101_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa'

# in_f='cnv.head.pvar'
# out_f='out.cnv.head.seq.pvar'
in_f=$1
out_f=$2
if [ $# -gt 2 ] ; then ref_fa=$3 ; fi

ml load R/3.6 gcc bedtools

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

# define intermediate file names
ucsc_bed_f="${tmp_dir}/tmp.$(basename ${in_f%.zst} .pvar).bed"
seq_f=${ucsc_bed_f%.bed}.tsv.gz

############################################################
## step 1: pvar --> UCSC bed format
############################################################

Rscript /dev/stdin ${in_f} ${ucsc_bed_f} << EOF
suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))

args <- commandArgs(trailingOnly=TRUE)
in_f <- args[1]
out_f <- args[2]

in_f %>% fread(colClasses=c('#CHROM'='character')) %>% rename('CHROM'='#CHROM') %>%
separate(ID, c('chr', 'POS_s', 'POS_e', NA), remove=F) %>%
mutate(indel = str_sub(ID, -1, -1)) %>%
select(chr, POS_s, POS_e, ID, indel) %>% rename('#chr' = 'chr') %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)
EOF

############################################################
## step 2: fetch sequence w/ bedtools getfasta
############################################################

bedtools getfasta -fi ${ref_fa} -bed ${ucsc_bed_f} -fo /dev/stdout -tab -name | bgzip > ${seq_f}

############################################################
## step 3: join the tables and write the pvar file
############################################################

Rscript /dev/stdin ${in_f} ${seq_f} ${out_f} << EOF
suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))

args <- commandArgs(trailingOnly=TRUE)
in_f <- args[1]
seq_f <- args[2]
out_f <- args[3]

seq_f %>% fread(col.names=c('ID', 'seq'), head=F) %>%
mutate(seq1 = str_sub(seq, 1, 1)) -> seq_df

in_f %>% fread(colClasses=c('#CHROM'='character')) %>% rename('CHROM'='#CHROM') %>%
select(-REF, -ALT, -POS) %>%
separate(ID, c(NA, 'POS', NA, NA), remove=F) %>%
mutate(indel = str_sub(ID, -1, -1)) %>% left_join(seq_df, by='ID') %>%
mutate(REF=if_else(indel == '+', seq1, seq), ALT=if_else(indel == '+', seq, seq1)) %>%
select(CHROM, POS, ID, REF, ALT) %>% 
arrange(CHROM, POS, ID) %>% rename('#CHROM' = 'CHROM') %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)
EOF

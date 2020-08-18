#!/bin/bash
set -beEuo pipefail

ml load R/3.6 gcc

for pop in 'white_british' 'non_british_white' 'african' 's_asian' 'e_asian' 'related' 'others' ; do
# pop="white_british"

evec_f="/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200313/pca_20200817_v1/ukb24983_${pop}.eigenvec"
evec_p="${evec_f}.PC1.PC2.png"

Rscript /dev/stdin ${evec_f} ${evec_p} << EOF
suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))
args <- commandArgs(trailingOnly=TRUE)

evec_f <- args[1]
evec_p <- args[2]

evec <- fread(evec_f, colClasses=c('#FID'='character', 'IID'='character')) %>%
rename('FID'='#FID')

p <- evec %>% 
ggplot(aes(x=PC1, y=PC2)) +
stat_density_2d(aes(fill = ..level..), geom = "polygon") +
labs(title = file.path(basename(dirname(evec_f)), basename(evec_f))) +
theme_bw()

ggsave(evec_p, p)
EOF

echo ${evec_p}
done

#!/bin/bash

zcat $(ls /scratch/groups/mrivas/users/ytanigaw/20200114_FinnGen_R2/summary_stats_hg19/* | grep -v unmapped) | cut -f1,2,3,4 | sort -u >finngen.bim

#sed '1iCHROM\tPOS\tREF\tALT\tMAF' >finngen.bim

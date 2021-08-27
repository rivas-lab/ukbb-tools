#!/bin/bash

gunzip -c mrp_rv_ma_exome_gbe.tsv.gz | awk -F'\t' 'BEGIN{OFS="\t"}($13 >= 5 && $13 != "nan" || $15 >= 5 && $15 != "nan" || $26 >= 5 && $26 != "nan" ){print $1,$2, $3,$4,$5, $6,$17, $13,$15,$26}' > ../19_mrpmm/mrp200kexomeresults.tsv
gunzip -c mrp_rv_ma_array_gbe.tsv.gz | awk -F'\t' 'BEGIN{OFS="\t"}($13 >= 5 && $13 != "nan" || $15 >= 5 && $15 != "nan" || $26 >= 5 && $26 != "nan" ){print $1,$2, $3,$4,$5, $6,$17, $13,$15,$26}' > ../19_mrpmm/mrparrayresults.tsv
gunzip -c mrp_rv_ma_exome_ultrarare_gbe.tsv.gz | awk -F'\t' 'BEGIN{OFS="\t"}($13 >= 5 && $13 != "nan" || $15 >= 5 && $15 != "nan" || $26 >= 5 && $26 != "nan" ){print $1,$2, $3,$4,$5, $6,$17, $13,$15,$26}' > ../19_mrpmm/mrp200kexomeultrarareresults.tsv

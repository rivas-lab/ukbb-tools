#!/bin/bash

touch field_table.tsv;

for file in $(ls /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/*/*/download/*.tab.columns | grep -v current); do 
    awk '{split(FILENAME, a, "/"); if (NR > 1) {print a[9]"\t"$3}}' $file | sort -u >> field_table.tsv;
done

sort -k1 field_table.tsv >>field_tmp && mv field_tmp field_table.tsv

join -1 1 -2 2 field_table.tsv <(cut -f1,2,7 ../tables/ukbb24983_tables.tsv | tail -n +2 | sort -u | sort -k2) -t $'\t' >../tables/field_table_basket_date.tsv

echo -e "TableID\tFieldID\tBasketID\tRelease_Date" | cat - ../tables/field_table_basket_date.tsv > temp && mv temp ../tables/field_table_basket_date.tsv

rm field_table.tsv

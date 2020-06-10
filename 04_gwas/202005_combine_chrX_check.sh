#!/bin/bash

# output
out_lst=202005_combine_chrX_check.$(date +%Y%m%d-%H%M%S).lst

# input params
data_dir="/oak/stanford/groups/mrivas/ukbb24983/cal/gwas/current"
prefix="ukb24983_v2_hg19."
suffices=(.array-combined.glm.logistic.hybrid.gz .array-combined.glm.linear.gz)
pops=(
african
e_asian
non_british_white
s_asian
white_british
)
#pop="white_british"

for pop in ${pops[@]} ; do
for suffix in ${suffices[@]} ; do
ls ${data_dir}/${pop}/*${suffix} |
sed -e "s%${data_dir}/${pop}/${prefix}%%g" |
sed -e "s%${suffix}%%g" |
sed -e "s%_X$%%g" |
sort | uniq -c |
awk -v OFS='\t' -v pop=${pop} '($1 == 2){print pop, $2}'
done
done > ${out_lst}

echo ${out_lst}

exit 0

ls /oak/stanford/groups/mrivas/ukbb24983/cal/gwas/current/white_british/*glm.logistic.hybrid.gz | sed -e "s%/oak/stanford/groups/mrivas/ukbb24983/cal/gwas/current/white_british/ukb24983_v2_hg19.%%g" | sed -e "s%.array-combined.glm.logistic.hybrid.gz%%g" | sed -e "s/_X$//g" | uniq -c | awk '$1 != 2'


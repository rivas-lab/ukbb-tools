#!/bin/bash
set -beEuo pipefail

pop=$1
columnar=$2

ml load python/3.6

data_d="/scratch/groups/mrivas/ukbb24983/exome/gwas/freeze/master_phe_20201002_exomeOQFE_20201110"

in_f="${data_d}/ukb24983_exomeOQFE.${pop}.p1e-3.tsv.gz"
out_f="${in_f%.tsv.gz}.${columnar}"

unzip_in_f=${in_f%.gz}

if [ ! -s ${unzip_in_f} ] ; then
    pigz -dc ${in_f} | sed -e 's/^#//g' > ${unzip_in_f}
fi

python3 11_tsv2columnar.py ${unzip_in_f} ${out_f} ${columnar}

exit 0
##################################
# usage:
bash 11_tsv2columnar.sh e_asian parquet
bash 11_tsv2columnar.sh e_asian feather

for pop in 'white_british' 'non_british_white' 'african' 's_asian' 'related' 'others' ; do
    sbatch -p mrivas --qos=high_p --nodes=1 --mem=12000 --cores=2 --time=7-0:00:00 \
        --job-name=columnar.${pop} --output=logs/columnar.${pop}.%A.out --error=logs/columnar.${pop}.%A.err \
        11_tsv2columnar.sh ${pop} parquet
done


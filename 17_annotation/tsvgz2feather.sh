#!/bin/bash
set -beEuo pipefail

annot_tsvgz=$1

ml load python/3.6

python3 tsvgz2feather.py ${annot_tsvgz} ${annot_tsvgz%.tsv.gz}.feather

exit 0
#######################
bash tsvgz2feather.sh \
    /oak/stanford/groups/mrivas/ukbb24983/array-combined/annotation/annotation_20201012/ukb24983_cal_hla_cnv.annot_20201023.tsv.gz

sbatch -p mrivas --nodes=1 --mem=300000 --cores=2 --time=6:00:00 \
    --job-name=tsvgz2feather --output=logs/tsvgz2feather.%A.out --error=logs/tsvgz2feather.%A.err \
    tsvgz2feather.sh \
    /oak/stanford/groups/mrivas/ukbb24983/exome/annotation/20201025_exome_oqfe_2020/ukb24983_exomeOQFE.annotation.tsv.gz


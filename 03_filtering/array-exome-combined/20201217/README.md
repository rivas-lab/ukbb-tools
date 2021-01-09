# array & exome combined dataset (hg19)

## Yosuke Tanigawa

In some applications, we want to take the union of the array (array-combined dataset) and exome (200k) dataset.

We use GRCh37/hg19 coordinate (use the exome datasets after applying liftUnder - the coverage of the chain file is relatively good in the coding region)

| geno_data_source | n          |
|------------------|-----------:|
| exome200k        | 17,660,693 |
| cal              |    802,498 |
| cnv              |    275,180 |
| hla              |        328 |
| Total            | 18,738,699 |

## data location

- `/oak/stanford/groups/mrivas/ukbb24983/array-exome-combined/pgen/ukb24983_cal_hla_cnv_exomeOQFE.pvar.{gz,zst}`
- `/scratch/groups/mrivas/ukbb24983/array-exome-combined/annotation/20201217/ukb24983_cal_hla_cnv_exomeOQFE.annot_20210108.tsv.gz`

## version history

- version `20210108`
  - there is no change in the genotype matrix
  - we used the latest variant annotation file (exome 200k has the QC flag in the annotation file) and generated the annotation file.
- version `20201217`
  - the initial attempt of the array-exome combined dataset
  - since we had not finalized the variant QC for the exome 200k dataset back then, we have all the variants

## analysis notebook and scripts

- Scripts and notebooks to generate the merged genotype matrix file
  - [`1a_intersection_pvar.ipynb`](1a_intersection_pvar.ipynb): identify non-overlapping set of variants and save it as a pvar file
    - The file is saved under `/scratch/groups/mrivas/ukbb24983/array-exome-combined/pgen/merge_list_pvar`
  - [`1b_sort_pvar.sh`](1b_sort_pvar.sh): sort the pvar file
  - [`1c_merge.sh`](1c_merge.sh): we use plink1.9 and merge the genotype dataset
  - [`1d_bed2pgen.sh`](1d_bed2pgen.sh): we convert plink1.9 BED dataset into plink2.0 pgen dataset
  - [`1e_add_back_source.ipynb`](1e_add_back_source.ipynb): we add one column to pvar file to keep track of the origin of the genotype data.
- [`2a_combine_annotation.R`](2a_combine_annotation.R): combine the variant annotation file.
  - [`2a_combine_annotation_draft.ipynb`](2a_combine_annotation_draft.ipynb): draft notebook
- [`combine_array_and_exome.sh`](combine_array_and_exome.sh): given a two files (one from array and the other from exome), we combine them together and write it to a specified output file. Internally, this calls [`combine_array_and_exome.R`](combine_array_and_exome.R).

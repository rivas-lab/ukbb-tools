# array & exome combined dataset (hg19)

## version 2020/12/16-17, Yosuke Tanigawa

In some applications, we want to take the union of the array (array-combined dataset) and exome (200k) dataset.

| geno_data_source | n          |
|------------------|-----------:|
| exome200k        | 17,660,693 |
| cal              |    805,426 |
| cnv              |    275,180 |
| hla              |        362 |
| Total            | 18,741,661 |

## data location

- `/oak/stanford/groups/mrivas/ukbb24983/array-exome-combined/pgen/ukb24983_cal_hla_cnv_exomeOQFE.pvar.{gz,zst}`
  - This file has non-overlapping set of variants

## analysis notebook and scripts

- [`1a_intersection_pvar.ipynb`](1a_intersection_pvar.ipynb): identify non-overlapping set of variants and save it as a pvar file
- [`1b_sort_pvar.sh`](1b_sort_pvar.sh): sort the pvar file
- [`combine_array_and_exome.sh`](combine_array_and_exome.sh): given a two files (one from array and the other from exome), we combine them together and write it to a specified output file. Internally, this calls [`combine_array_and_exome.R`](combine_array_and_exome.R).
- [`2a_combine_annotation.R`](2a_combine_annotation.R): combine the variant annotation file.
  - [`2a_combine_annotation_draft.ipynb`](2a_combine_annotation_draft.ipynb): draft notebook

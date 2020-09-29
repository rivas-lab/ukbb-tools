## merge array + imputation for LD map computation

We merge array + HLA + impuation for the [LD map computation](/14_LD_map/array_imp_combined_no_cnv).

Note: The CNVs were not included in this merge.

## File location

- `/oak/stanford/groups/mrivas/ukbb/24983/array_imp_combined_no_cnv/pgen`
- We have two versions of the file. Please see the version history section below.
  - `pgen_v1`
  - `pgen_v2`


## Version history

## `pgen_v2`

As we describe in [the variant QC for the imputed dataset](/03_filtering/imp/3_merge), we changed the variant QC criteria and decided NOT to use the HWE filter.

Our source datasets are listed in [`merge.lst.v2.imp.wo.HWE.filter.tsv`](merge.lst.v2.imp.wo.HWE.filter.tsv).

We subsequently generated the `array_imp_combined_no_cnv` dataset using the updated QC criteria as version 2 of the dataset.

## `pgen_v1`, 2019/11/21

We used version 3 of the imputed dataset and applied our variant QC criteria (version 1). 

Our source datasets are listed in [`merge.lst.v1.imp.w.HWE.filter.tsv`](merge.lst.v1.imp.w.HWE.filter.tsv).

We generated the `array_imp_combined_no_cnv` dataset. Please see [the variant QC for the imputed dataset](/03_filtering/imp/3_merge) for more information.


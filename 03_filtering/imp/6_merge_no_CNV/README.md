# merge array + imputation for LD map computation (no longer supported)

We used to have `cal` + `hla` + `imp` (QC-passed) genotype dataset for the [LD map computation](/14_LD_map/array_imp_combined_no_cnv).

We later realized that this dataset was redundant, because we could run the LD map comuptation using `cal` + `hla` + `cnv` + `imp` (QC-passed) genotype dataset by specifying the list of variants we'd like to include in the LD map.

## File location

- We used to have the genetic data here: `/oak/stanford/groups/mrivas/ukbb/24983/array_imp_combined_no_cnv/pgen`
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

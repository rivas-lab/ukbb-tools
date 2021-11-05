# Imputation v3 dataset variant QC

We apply the following variant QC procedure.

## Original variant QC procedure (version 1, no longer used)

The steps up to variant position comparison is documented in notebook (`imp_v3_QC.ipynb`) and HWE-based filter and extraction of those variants are in a script (`imp_extract.sh`).

```
97,059,329 Variants on the imputation dataset
  |
  | MAF >= 0.01
  | Imputation quality >= 0.7
  |  
  | where, MAF and Imputation quality is obtained from UKB data showcase
  | (http://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=1967)
  |
  | Note, MAF is computed for the entire 500k (by UKB)
  |
10,061,256
  |
  | Biallelic only
  |
10,028,119
  |
  | The variant position is not present on genotyping array dataset
  |
 9,354,600
  |
  | Missingness
  | HWE 1e-7
  |
 5,182,706
```

The list of 9,354,600 variants is recorded in: `/oak/stanford/groups/mrivas/ukbb24983/imp/mfi/ukb_mfi_v3.info.maf.biallelic.noncal.tsv.zst`

## Alternative stretegy without HWE filter (version 2)

Similar, we performed a very similar filtering except for the HWE filter at the end. [`imp_extract_noHWE.sh`](imp_extract_noHWE.sh) was used for this analysis. The resulting pgen files (per chromosome) are in `/oak/stanford/groups/mrivas/ukbb/24983/imp/pgen/var_QC_noHWE/ukb24983_imp_chr${chr}_v3_QC_noHWE.[pgen,psam,pvar.zst]`.

```
97,059,329 Variants on the imputation dataset
  |
  | MAF >= 0.01
  | Imputation quality >= 0.7
  |  
  | where, MAF and Imputation quality is obtained from UKB data showcase
  | (http://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=1967)
  |
  | Note, MAF is computed for the entire 500k (by UKB)
  |
10,061,256
  |
  | Biallelic only
  |
10,028,119
  |
  | The variant position is not present on genotyping array dataset
  |
 9,354,600
  |
  | Missingness
  |
 5,931,362
```

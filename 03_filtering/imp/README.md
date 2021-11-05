# imputation dataset (imp) version 3

### Yosuke Tanigawa (`yosuke.tanigawa@alumni.stanford.edu`)

This directory contains some pre-processing scripts for the imputation dataset (version 3).

Briefly, we
1. converted the BGEN file from UK Biobank into plink2 pgen format;
2. Variant QC (please see the `2_var_QC` dir for more details):
3. generated "array_imp_combined" merge file;
4. generated training/validation/test dataset;
5. apply flip_check for the resulting dataset;

## Two versions of the variant QC criteria

In the initial variant QC (Dec. 2019), we used the HWE filter.

In May 2020, we revised this filter and decided NOT to use this HWE filter.

Please read more in [`2_var_QC`](2_var_QC). The 

## File location

- `/oak/stanford/groups/mrivas/ukbb24983/imp/pgen`
  - This directory should have the imputed genotype data (before variant QC) in the plink 2 format.
- `/oak/stanford/groups/mrivas/ukbb24983/array_imp_combined/pgen_v2`
  - The genotype data in this directory has the combination of the followigns: the directly genotyped variants, imputed HLA allelotypes, CNVs, and the QC-passed imputed variants that do not overlap with the directly genotyped variants.
- `/oak/stanford/groups/mrivas/ukbb24983/imp/annotation`
  - Variant annotation from Ensembl's VEP for the imputed dataset

## References

- https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=263

# imputation dataset (imp) version 3

This directory contains some pre-processing scripts for the imputation dataset (version 3).

Briefly, we
1. converted the BGEN file from UK Biobank into plink2 pgen format;
2. Variant QC (please see the `2_var_QC` dir for more details):
3. generated "array_imp_combined" merge file;
4. generated training/validation/test dataset;
5. apply flip_check for the resulting dataset;

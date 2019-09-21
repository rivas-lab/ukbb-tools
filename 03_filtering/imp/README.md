# imputation dataset (imp) versino 3

This directory contains some pre-processing scripts for the imputation dataset (version 3).

Briefly, we
1. converted the BGEN file from UK Biobank into plink2 pgen format;
2. identified biallelic variants with MAF >= 1 % ;
3. generated "array_imp_combined" merge file;
4. generated training/validation/test dataset.

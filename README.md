# ukbb-tools

This repository contains our complete set of tools for preprocessing, quality control, and preliminary analyses on UK Biobank data. There is a folder in the repo per set of methods as defined in the Table of Contents below. Each subdirectory has a README.md file that should be read before use. These files detail how to use all files within the directory.

## Contents
1. [Preprocessing](https://github.com/rivas-lab/ukbb-tools/tree/master/01_preprocessing)
2. [Phenotyping](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping)
3. [Filtering](https://github.com/rivas-lab/ukbb-tools/tree/master/03_filtering)
4. [GWAS](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas)
5. [GBE](https://github.com/rivas-lab/ukbb-tools/tree/master/05_gbe)
6. [PheWAS](https://github.com/rivas-lab/ukbb-tools/tree/master/06_phewas)
7. [LD score regression (LDSC)](https://github.com/rivas-lab/ukbb-tools/tree/master/07_LDSC)
8. [UK Biobank Bulk Download](https://github.com/rivas-lab/ukbb-tools/tree/master/08_bulk_DL)
9. [Coordinate Liftover (UCSC)](https://github.com/rivas-lab/ukbb-tools/tree/master/09_liftOver)
10. [Biomarker Adjustment](https://github.com/rivas-lab/ukbb-tools/tree/master/10_phe_adjustment)
11. [GREAT Enrichment Calculator](https://github.com/rivas-lab/ukbb-tools/tree/master/11_great_enrichment)
12. [SciDB Query for PheWAS](https://github.com/rivas-lab/ukbb-tools/tree/master/12_query_scidb)
13. [Multiple Rare-variants and Phenotypes (MRP) - Rare-variant signal aggregator](https://github.com/rivas-lab/ukbb-tools/tree/master/13_mrp)
14. [LD map](https://github.com/rivas-lab/ukbb-tools/tree/master/14_LD_map)
15. [Genetic Relationship Matrix calculation (GRM via GCTA)](https://github.com/rivas-lab/ukbb-tools/tree/master/15_GRM)
16. [snpnet (Large-scale Cox Proportional Hazards)](https://github.com/rivas-lab/ukbb-tools/tree/master/16_snpnetcox)
17. [VEP Variant Annotation](https://github.com/rivas-lab/ukbb-tools/tree/master/17_annotation)
18. [METAL (Meta-analysis)](https://github.com/rivas-lab/ukbb-tools/tree/master/18_metal)
19. [Multiple Rare-variants and Phenotypes Mixed Model (MRPMM)](https://github.com/rivas-lab/ukbb-tools/tree/master/19_mrpmm)

## The `ukbb-tools` module on Sherlock

All this code has been ported to a module on Sherlock. [Click](https://github.com/rivas-lab/sherlock-modules/tree/master/ukbb-tools) for more details on how to load and use this module.

There is an updater script that pushes your current directory - use with appropriate caution, as it takes the `master` branch - and makes it a version of the module. The only argument for the updater is a date; this is used as a version label.

Example Usage:
```{bash}
$ bash ukbb-tools.module.updater.sh 20200225
``

## Allele frequency across UKB populations in the array-combined dataset

`/oak/stanford/groups/mrivas/ukbb24983/array-combined/afreq_20201012/ukb24983_cal_hla_cnv.afreq_20201012.pvar.zst`

This (zstd-compressed) table file has 73 columns.

- The first 5 columns represent the location and the alleles of the variant as in the `pvar` file
  - `CHROM`, `POS`, `ID`, `REF`, and `ALT`
- Then, we have `array` column representing whether the variant is genotyped on BiLEVE array (`UKBL`), Axiom array (`UKBB`), or both (`both`). For CNVs and HLA allelotypes, we place `NA`.
- The next 64 columns are allele frequency and allele count computed for (1 + 7) populations in UK Biobank.
  - The (1 + 7) populations are:
    - The full 488k dataset of genotyped individuals
    - The 7 populations as defined in [`population_stratification_20200828`](03_filtering/population_stratification_20200828)
      - `white_british`, `non_british_white`, `african`, `s_asian`, `e_asian`, `related`, `others`
  - The 8 characteristics are:
    - `AF`: Alternate allele frequency, computed from `plink2 --freq` command
    - `OBS_CT`: the number of individuals with non-missing genotype values (`plink2 --geno-counts cols=nobs`)
    - `MISSING_CT`: the number of individuals with missing genotype values (`plink2 --geno-counts cols=missing`)
    - `HOM_REF_CT`: the number of individuals with homozygous reference allele (`plink2 --geno-counts cols=homref`)
    - `HET_REF_ALT_CTS`: the number of individuals with heterozygous allele (`plink2 --geno-counts cols=refalt`)
    - `TWO_ALT_GENO_CTS`: the number of individuals two alternate alleles (this dataset should represent biallelic variabts, so this should be the number of individuals with homozygous alternate allele) (`plink2 --geno-counts cols=altxy`)
    - `HAP_REF_CT`: the number of individuals with hemizygous reference allele (`plink2 --geno-counts cols=hapref`)
    - `HAP_ALT_CTS`: the number of individuals with hemizygous alternate allele (`plink2 --geno-counts cols=hapalt`)
- The next 3 columns represent the missingness 
    - `f_miss`: missing rate
    - `f_miss_UKBL`: missing rate in UKBL array
    - `f_miss_UKBB`: missing rate in UKBB array

### Note

- we used the latest master sqc file to get the number of individuals genotyped on each array.
  - master sqc file `/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200828/ukb24983_master_sqc.20200828.phe`
  - {`UKBB`: 438427, `UKBL`: 49950}

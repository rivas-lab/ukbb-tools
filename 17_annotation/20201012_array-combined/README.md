# variant annotation and variant QC for the array-combined dataset

### Yosuke Tanigawa (work in progress)

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


## HWE

For HWE test in chrX, plink implements this relatively new procedure.

- J. Graffelman, B. S. Weir, Testing for Hardy–Weinberg equilibrium at biallelic genetic markers on the X chromosome. Heredity. 116, 558–568 (2016). https://doi.org/10.1038/hdy.2016.20

This method does NOT ignore the males when applying HWE test for chr X.

And here is the HWE p-value distribution (the red bar indicates our P-value cutoff of 1e-7):

![HWE midp plot](hwe_midp_plot.png)

## Variant QC

Here is the summary of variant QC (across both autosomal and non-autosomal variants in the combined array dataset).

![variant QC summary](variant_QC.png)

We have the following QC filters: 

- `missingness`: missingness (1%, computed separately for UKBL/UKBB array if the variant is directly genotyped and present in only one array)
- `hwe`: HWE p-value (1e-7, we used the chrX model above)
- `mcpi`: manual cluster plot inspection (copied from variant QC file back in 2017)
- `gnomad_af`: maf comparison with gnomAD  (copied from variant QC file back in 2017)
- `mgi`: manual genome browser inspection???  (copied from variant QC file back in 2017)


## LD pruning

We apply LD pruning with `plink2 --indep-pairwise 1000kb 1 0.5`.
Note, we were originally using `plink2 --indep-pairwise 50 5 0.5`, but we swapped to distance-based LD window specification so that we are consistent with what we do in the LD map comptuation.

To prioritize the variants with severe predicted consequence, we applied LD pruning for each consequence group using the following procedure.

0. We focus on the QC-passed genotyped variants for this LD pruning analysis. In other words, we don't use CNVs and HLA allelotypes for the LD pruning analysis.
1. We apply LD pruning for the protein-truncating variants (PTVs). We check the pre-computed LD map with the same `r2` threshold and remove the variants that are in linkage with the selected variants.
2. We apply LD pruning for the ptotein-altering variants (PAVs) that are not in LD with previously selected variants (PTVs). Using the union of the selected PTVs and PAVs, we check the LD map and remove the variants that are in linkage.
3. Similarly, we apply LD pruning for the protein-coding variants (PCVs) that are no in LD with the previously selected PTVs or PAVs. Using the union of the selected PTVs, PAVs, and PCVs, we check the LD map and remove the variants in linkage.
4. Repeat the procedure for UTR-region variants, intronic variants, and the remaining variants.




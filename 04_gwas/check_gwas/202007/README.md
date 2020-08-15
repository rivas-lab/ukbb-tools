# checking GWAS output files (2020 spring/summer freeze)

In this freeze, we created GWAS QC table.

The results are in `/oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current/gwas.qc.tsv` and we have [a copy in Google spreadsheet](https://docs.google.com/spreadsheets/d/12cQ9NQcj5jRWY4VZtxMU9E01drUc_grUuNbV_UGiD8o/edit?usp=sharing).

### column descriptor:

The table file has the following set of columns.

- Phenotype and population names:
  - `GBE_ID`: the GBE phenotype code
  - `population`: the analyzed population. `metal` represents the UKB-wide meta-analysis
- Phenotype information
  - `GBE_category`: the phenotype category shown in GBE
  - `N`: the number of cases (binary traits) or the number of non-NA individuals (quantitative traits) used in the GWAS analysis
  - `GBE_NAME`: the full phenotype name
  - `GBE_short_name`: the short phenotype name
- The number of lines and the number of genome-wide significant (p < 5e-8) hits
  - `n_lines`
  - `n_non_NA_lines`
  - `n_hits`
  - `n_ld_indep_hits`
- The lambda GC statistics across the frequency bins
  - `lgc.0001`
  - `lgc.0001-.001`
  - `lgc.001-.01`
  - `lgc.01`
  - `lgc.05`
  - `lgc.common`
- The statistics from LDSC heritability estimate
  - `h2_obs`
  - `h2_obs_se`
  - `lambda_GC`
  - `mean_chi2`
  - `intercept`
  - `intercept_se`
  - `ratio`
  - `ratio_se`

## Methods

Issue ticket corresponding to this project: https://github.com/rivas-lab/ukbb-tools/issues/32

### Phenotype information

We took the phenotype category information from [`05_gbe`](/05_gbe/extras/20200812_GBE_category).

### Line counts and the number of hits

- [`line.count.sh`](line.count.sh): this script counts the number of lines and the number of hits. Specifically, this script takes the GWAS summary statistics (in `/gwas/current/<population>` directory) and generate a small table consists of the following columns:
  - `GBE_ID`: The phenotype identifier we use in the lab (and also on GBE)
  - `population`: The population name string. The `metal` represents the UKB-wide meta-analysis results
  - `n_lines`: The number of lines in the summary statistics file
  - `n_non_NA_lines`: The number of non-NA (based on the P-value column) lines in the summary statistics file
  - `n_hits`: The number of genome-wide significant hits (univariate P-value < 5e-8) in the summary statistics file
  - `n_ld_indep_hits`: The number of LD independent genome-wide significant hits (univariate P-value < 5e-8) in the summary statistics file. The list of LD independent variants are taken from the variant annotation file for the array dataset.

We compute those statistics for each summary statistics file and aggregate them into one table (per population) using [`line.count.agg.sh`](line.count.agg.sh).

### lambda GC across different frequency bins 

- [`computelqc.R`](computelqc.R): this script computes the Lambda GC across the following frequency bins

| frequency bin string | the maf filter             |
|----------------------|----------------------------|
| common               | maf > .05                  |
| 0.05                 | maf <= .05                 |
| 0.01                 | maf <= .01                 |
| .001-.01             | maf <= .01 & maf >= .001   |
| .0001-.001           | maf <= .001 & maf >= .0001 |
| 0.0001               | maf <= .0001               |

We compute those statistics for each summary statistics file and aggregate them into one table (per population) using [`computelqc.agg.sh`](computelqc.agg.sh).

### LDSC intercept and heritability estimates

We took the results from the recent LDSC h2 analysis. Briefly, we applied `ldsc.py --h2` using our [wrapper script](07_LDSC/helpers/ldsc_h2.sh). Please check our [LDSC job scripts](/07_LDSC/jobs/202007_LDSC) and [LDSC official documentation](https://github.com/bulik/ldsc) for more information.

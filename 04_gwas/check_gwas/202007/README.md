# checking GWAS output files (2020 spring/summer freeze)

In this freeze, we created GWAS QC table containing the following metrics:

- Lambda GC across different frequency bins
- LD score regression (LDSC)
  - heritability estimate
  - LDSC intercept
- Line counts
  - Number of lines
  - Number of non-NA lines
- The number of hits
  - Number of hits (p < 5e-8)
  - Number of independent hits (load in LD pruned set)

Issue ticket corresponding to this project: https://github.com/rivas-lab/ukbb-tools/issues/32

## Methods

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

### Line counts and the number of hits

- [`line.count.sh`](line.count.sh): this script counts the number of lines and the number of hits. Specifically, this script takes the GWAS summary statistics (in `/gwas/current/<population>` directory) and generate a small table consists of the following columns:
  - `GBE_ID`: The phenotype identifier we use in the lab (and also on GBE)
  - `population`: The population name string. The `metal` represents the UKB-wide meta-analysis results
  - `n_lines`: The number of lines in the summary statistics file
  - `n_non_NA_lines`: The number of non-NA (based on the P-value column) lines in the summary statistics file
  - `n_hits`: The number of genome-wide significant hits (univariate P-value < 5e-8) in the summary statistics file
  - `n_ld_indep_hits`: The number of LD independent genome-wide significant hits (univariate P-value < 5e-8) in the summary statistics file. The list of LD independent variants are taken from the variant annotation file for the array dataset.

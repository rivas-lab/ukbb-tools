# extras - miscellaneous scripts related to GWAS analysis

## helper scripts

### [`gwas_filter_by_p.R`](gwas_filter_by_p.R)

#### Author: Yosuke Tanigawa

We often have many rows in the GWAS summary statistics and are generally interested in the subset of them. We typically use the P-value threshold and the summary statistics with `P > 1e-3`, for example are not used in most of the analyses.

With those observation in mind, we prepare summary statistics file with p-value filter. One naive way is to use a command line tool like `awk`, but we learned that that may not work great, when we have `P << 1e-300` (below machine epsilon in double precision).

For that purpose, we have an R script called [`gwas_filter_by_p.R`](gwas_filter_by_p.R). This reads P-value column as a string and apply p-value threshold on log10(p) value.

#### usage

```
Rscript gwas_filter_by_p.R input.sumstats.tsv output.sumstats.tsv [p_threshold]
```

The `p_threshold` is an optional argument. If it's not specified, the script uses its default threshold of 1e-3.

### [`sumstats_to_plink.sh`](sumstats_to_plink.sh)

#### Author: Yosuke Tanigawa

This script converts a given summary statistics into plink format.

The usage is explained in its help message.

```{bash}
bash sumstats_to_plink.sh
sumstats_to_plink.sh: 1 positional arguments are required
sumstats_to_plink.sh (version 0.0.1)
Convert a summary statistics into a PLINK format.

Usage: sumstats_to_plink.sh [options] input_file
  input_file      The input file

Options:
  --logit   (-l)  Specify it as a logistic regression.

Default configurations:
  logit=FALSE
  col_key_CHROM=CHROM
  col_key_POS=POS
  col_key_ID=ID
  col_key_REF=REF
  col_key_ALT=ALT
  col_key_A1=A1
  col_key_OBS_CT=OBS_CT
  col_key_BETA=BETA
  col_key_OR=OR
  col_key_SE=SE
  col_key_P=P
```

## A primer on formatting external summary stats for GBE

This directory serves as a catch-all for additional summary statistics which should be added to the analysis pipeline and for inclusion into the [Global Biobank Engine](biobankengine.stanford.edu). This specification is general enough to work with any input source, but this documentation will operate under the assumption that the summary stats are derived from UK Biobank resources.

*Important*: This directory should _not_ be used for data where we have phenotype-level data from UK Biobank. If we have the raw phenotype files, you should follow the steps outlined in [`02_phenotyping/extras`](../../02_phenotyping/extras) instead.

### Specification

Each new data source should be tracked with its own subdirectory, kept in this folder. That folder should contain the following:

1. A script which processes raw summary statistics: This script will need to map and subset variants in the input to array-genotyped variants in the UK Biobank, and output the resulting files according PLINK's [specification](https://www.cog-genomics.org/plink/2.0/formats#glm_logistic) (bolded fields only is fine). This is currently (April 2019) a requirement for upload into GBE, as only array-genotyped sites present in UK Biobank function with the data loader.

2. A tab-delimited info file, called `info.tsv`: This must contain phenotype IDs in the first column, N in the second column, phenotype names in the third column, and paths to output files (generated in step 1) in the fourth column. See below for an example.

### Current extra GWAS results we have

1. [`MVP`](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/extras/MVP)
2. [`broad_collaboration`](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/extras/broad_collaboration)
3. [`bbj`](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/extras/bbj)

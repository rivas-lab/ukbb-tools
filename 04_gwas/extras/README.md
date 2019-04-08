# A primer on formatting external summary stats for GBE

This directory serves as a catch-all for additional summary statistics which should be added to the analysis pipeline and for inclusion into the [Global Biobank Engine](biobankengine.stanford.edu). This specification is general enough to work with any input source, but this documentation will operate under the assumption that the summary stats are derived from UK Biobank resources.

*Important*: This directory should _not_ be used for data where we have phenotype-level data from UK Biobank. If we have the raw phenotype files, you should follow the steps outlined in [02_phenotyping/extras](../../02_phenotyping/extras/README.md) instead.

## Specification
Each new data source should be tracked with its own subdirectory, kept in this folder. That folder should contain the following:

1. A script which processes raw summary statistics: This script will need to map and subset variants in the input to array-genotyped variants in the UK Biobank, and output the resulting files according PLINK's [specification](https://www.cog-genomics.org/plink/2.0/formats#glm_logistic) (bolded fields only is fine). This is currently (April 2019) a requirement for upload into GBE, as only array-genotyped sites present in UK Biobank function with the data loader.

2. A tab-delimited info file, called `info.tsv`: This must contain phenotype IDs in the first column, N in the second column, phenotype names in the third column, and paths to output files (generated in step 1) in the fourth column. See below for an example.

## Example: Summary Statistics from Neale Lab (Broad)
A readme on the project, as described by the Neale Lab can be found [here](http://www.nealelab.is/uk-biobank).

Numbers below correspond to the above requirements, followed by snippets of the actual files:
1. [project_to_biobank.py](broad_collaboration/project_to_biobank.py): 

2. [info.tsv](broad_collaboration/info.tsv):

```preview goes here```

Note that the result files are kept in `$OAK/ukbb24983/cal/gwas/extras/broad_collaboration`. Raw downloads are in `$OAK/ukbb/broad_collaboration`.

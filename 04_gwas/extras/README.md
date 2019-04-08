# A primer on formatting external summary stats for GBE

This directory serves as a catch-all for additional summary statistics which should be added to the analysis pipeline and for inclusion into the [Global Biobank Engine](biobankengine.stanford.edu). This specification is general enough to work with any input source, but this documentation will operate under the assumption that the summary stats are derived from UK Biobank resources.

*Important*: This directory should _not_ be used for data where we have phenotype-level data from UK Biobank. If we have the raw phenotype files, you should follow the steps outlined in [02_phenotyping/extras](../../02_phenotyping/extras/README.md) instead.

## Specification
Each new data source should be tracked with its own subdirectory, kept in this folder. That folder should contain the following:

1. A tab-delimited list of summary statistics, containing their location (paths on Sherlock), case counts (N), and name/ID for the Global Biobank Engine. Note that the raw data should be kept in the appropriate folder (see below) in the lab `$OAK` space.

2. A script which maps variants in the input summary statistics to overlapping variants in the UK Biobank. This is currently (April 2019) a requirement for upload into GBE, as only array-genotyped sites present in UK Biobank function with the data loader.

3. A minimal version of (1), called `info.tsv`, containing GBE IDs in the first column, N in the second column, and GBE names in the third column. See below for an example.

## Example: Summary Statistics from Neale Lab (Broad)
A readme on the project, as described by the Neale Lab can be found [here](http://www.nealelab.is/uk-biobank).

Numbers below correspond to the above requirements, followed by snippets of the actual files:
1. [neale_info.tsv](broad_collaboration/neale_info.tsv): 

Note that these files are kept in `$OAK/ukbb24983/imp/gwas/extras/broad_collaboration`.

2. [project_to_biobank.py](broad_collaboration/project_to_biobank.py): 

3. [info.tsv](broad_collaboration/info.tsv):

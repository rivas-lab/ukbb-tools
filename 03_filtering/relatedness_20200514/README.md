# The relatedness filter analysis

Hakhamanesh provided an initial analysis with the King's relationship table. His analysis notebook is included in this repo: [`UKB_relatedness.Rmd`](UKB_relatedness.Rmd). Sadly, GitHub does not render the R markdown document, so we have a redundant pdf copy of it: [`UKB_relatedness.pdf`](UKB_relatedness.pdf).

Note that he applied his own QC on individuals (not in this code) filtering individuals with high missingness, gender not matching genetically inferred sex, etc. After performing the GWAS power analysis, we will feed it your own standard set of QC’d individuals for consistency across analyses.

For now, we have the list of individuals in `/oak/stanford/groups/mrivas/ukbb24983/sqc/relatedness_20200514`.
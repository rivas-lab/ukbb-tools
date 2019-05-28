This directory contains documentation and scripts for additional phenotype definitions (sessions) that were processed in an ad-hoc fashion because of the nature of the data and/or the way they were packaged by UK Biobank. These projects include:

- Additional medication information
- Biomarkers from Basket 2001440, adjusted for covariates and/or statins
- Cancer diagnoses
- Clinical breathing Indices
- Phenotypes indicative of family history
- High-confidence phenotypes verified by multiple criteria
- Physical activity phenotypes obtained from Ashley lab
- Rohit's phenotypes

Each one of the above has its own script in its directory that makes the `.info` files for those phenotypes (and a symlink to `../../scripts/annotate_phe.py` in order to import the `make_phe_info` function). Each directory `XXXX` also has a `XXXX_gbe_map.tsv` that maps the names of the phenotypes to the GBE IDs.

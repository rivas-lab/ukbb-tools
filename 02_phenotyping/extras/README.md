# Phenotyping extras

## Contents

This directory contains subdirectories containing documentation and scripts for additional phenotypess that were processed in an *ad hoc* fashion because of the nature of the data and/or the way they were packaged by UK Biobank. These phenotypes are not readily processed by our [pipeline](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping#generating-and-updating-phenotypes-and-summary-statistics)

These projects include:

1. [Additional medication information](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/extras/additional_medications)
2. [Biomarkers from Basket 2001440, adjusted for covariates and/or statins](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/extras/adjusted_biomarkers)
3. [Cancer diagnoses](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/extras/cancer)
4. [Clinical breathing indices](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/extras/clinical_breathing_indices)
5. [Country of birth](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/extras/country_of_birth)
6. [Phenotypes indicative of family history](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/extras/family_history)
7. [High-confidence phenotypes verified by multiple criteria (and/or generated from first-occurrence/algorithmically defined criteria)](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/extras/highconfidenceqc)
8. [Additional intra-ocular_pressure phenotypes](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/extras/iop) for the [glaucoma project](https://github.com/rivas-lab/ANGPTL7/)
9. [Physical activity phenotypes obtained from Ashley lab](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/extras/physical_activity)
10. [Rohit's phenotypes](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/extras/rohit)
11. [Self-reported ethnicity phenotypes](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/extras/self_reported_ethnicity)

Each one of the above directories has its own `.py` script in its directory that makes the `.info` files for those phenotypes (and a symlink to [`../../scripts/annotate_phe.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/annotate_phe.py) in order to import the `make_phe_info` function). Each directory `XXXX` also has a `XXXX_gbe_map.tsv` that maps the names of the phenotypes to their GBE IDs, which have been assigned specifically to not conflict with existing phenotypes (see the current [allocation of phenotype names here](https://github.com/rivas-lab/wiki/blob/master/ukbb/icdinfo_allocation.md)).

## Pipelines and Workflows

Once phenotype and info files for your project are generated, and it is determined that these should then be routinely analyzed in the larger [pipeline](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping#generating-and-updating-phenotypes-and-summary-statistics), the [`symlink_extras.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/extras/symlink_extras.py) should be edited to include the source directory of your phenotypes (e.g. after line 51 in [`symlink_extras.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/extras/symlink_extras.py), add `your_dir = "/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/XXXX/"`, and edit line 53 to include `your_dir` as well), and step 9 of the [pipeline](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping#generating-and-updating-phenotypes-and-summary-statistics) should be run so that the relevant reference files are updated and GWAS can be run on your new phenotypes.

Please also look at the [`ukbb-phenotyping`](https://github.com/rivas-lab/ukbb-phenotyping) repository for the original sources of HC, death registry, and ICD/MED-related activities.

## SQL

Starting 2020, UKB also released some of the phenotypic data via SQL portal. This includes the covid-19 test results as well as the linked death registry information.

For instructions on to how to access this data see: http://biobank.ndph.ox.ac.uk/showcase/showcase/docs/DeathLinkage.pdf

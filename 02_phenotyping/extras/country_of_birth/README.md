# Country of birth information from UK Biobank

Yosuke Tanigawa

2020/3/28

We used the following two fields to extract the information on the country of birth.

- [Data-Field 20115: Country of Birth (non-UK origin)](http://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20115) ([Data-coding 89](http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=89))
- [Data-Field 1647: Country of birth (UK/elsewhere)](http://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1647) ([Data-coding 100420](http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100420))

We focued on individuals with array genotype data, removed the individuals with participation withdrawal, and applied the following QC filters in sample QC file.

1. `putative_sex_chromosome_aneuploidy` (n = 652)
2. `het_missing_outliers` (n = 968)

We additionally recoreded whether the individual is in the maximal set of unrelated individuals, which is identified by UK Biobank (`pass_filter` column).

We investegated the answers for the country of birth question and focused on individuals with consistent answer after removing uninformative options (such as do not know and prefer not to answer).

We extracted the following three types of information:

1. `country_of_birth_id` and `country_of_birth_txt`: this one is a combination of the following two information, (2) and (3). Specifically, we replaced the raw finest grade information with the corresponding top-level one when the total number of individuals is less than 30 within the set of 488k individuals.
2. `country_of_birth_top_level_id` and `country_of_birth_top_level_txt`: this one records the top-level classification of the country of birth information.
3. `country_of_birth_raw_id`  and `country_of_birth_raw_txt`: this one is records the finest information available from the fields.

Note that the id string is a contatenation of the field ID and the data-coding string with `_` as separater. Please see [`1_country_of_birth_QC.ipynb`](1_country_of_birth_QC.ipynb) for more information.

## data location

- Extracted file
  - `/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/country_of_birth/misc/ukb2005693_ukb37855.fields1647_20115.tsv.zst`
- Data-Coding mapping files
  - `coding89.tsv`
  - `coding100420.tsv`
- Counts of individuals per country/region
  - `/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/country_of_birth/misc/ukb2005693_ukb37855.country_of_birth.coding.tsv`

## scripts and notebook

- [`extract_country_of_birth.R`](extract_country_of_birth.R): this R script extracts the relevant data field from a tab file.
- [`1_country_of_birth_QC.ipynb`](1_country_of_birth_QC.ipynb): this notebook applies QC filters.

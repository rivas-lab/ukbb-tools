# array-combined dataset

- `both_array_list.sh`: the variants in `cal` dataset is measured in two genotyping arrays.


## penalty factor for snpnet

We have penalty factor files for snpnet.

The version 1 is pepared by Manny, version 2-5 are prepared by Yosuke.

### version 5

  - `/oak/stanford/groups/mrivas/ukbb24983/array-combined/snpnet/penalty.v5.rds`
  - New VEP run (version 101, 2020 Aug). We curated its consequence field and grouped into categories, such as PTVs and PAVs. Please read [this page](/17_annotation/20200912_cal_vep_v101).
  - We considered ClinVar information.
    - [`clinvar_extract.sh`](clinvar_extract.sh): extract the pathogenicity info from ClinVar vcf file
    - [`array-combined.p_factor.v5.R`](array-combined.p_factor.v5.R): clean-up ClinVar pathogenicity string and generate penalty factor file.
  - p.factor assignment rule
    - PTVs or ClinVar pathogenic --> 0.5
    - PAVs or ClinVar Likely pathogenic --> 0.75
    - Others --> 1.0
 
| p.factor | consequence | ClinVar           | n      |
|----------|-------------|-------------------|--------|
| 0.5      | ptv         |                   | 23494  |
| 0.5      | ptv         | Pathogenic        | 4421   |
| 0.5      | pav         | Pathogenic        | 3282   |
| 0.5      | ptv         | Likely_pathogenic | 406    |
| 0.5      | intron      | Pathogenic        | 51     |
| 0.5      | pcv         | Pathogenic        | 12     |
| 0.5      | others      | Pathogenic        | 7      |
| 0.75     | pav         |                   | 84937  |
| 0.75     | pav         | Likely_pathogenic | 942    |
| 0.75     | pcv         | Likely_pathogenic | 11     |
| 0.75     | intron      | Likely_pathogenic | 10     |
| 0.75     | others      | Likely_pathogenic | 1      |
| 1        | intron      |                   | 358378 |
| 1        | others      |                   | 310287 |
| 1        | pcv         |                   | 11259  |
| 1        | utr         |                   | 7928   |

### version 4
  - `/oak/stanford/groups/mrivas/ukbb24983/array-combined/snpnet/penalty.v4.rds`
  - New VEP run (version 101, 2020 Aug). We curated its consequence field and grouped into categories, such as PTVs and PAVs. Please read [this page](/17_annotation/20200912_cal_vep_v101).
  - We found that having too many variants with `w=.9` were not helpful.

| Data source | Csq    | n      | Csq_priority | w    |
|-------------|--------|--------|--------------|------|
| Array       | ptv    | 28321  | 1            | 0.5  |
| Array       | pav    | 89161  | 2            | 0.75 |
| HLA         |        | 362    |              | 0.75 |
| Array       | pcv    | 11282  | 3            | 1    |
| Array       | intron | 358439 | 4            | 1    |
| Array       | utr    | 7928   | 5            | 1    |
| Array       | others | 310295 | 6            | 1    |
| CNVs        |        | 275180 |              | 1    |

### version 3
  - `/oak/stanford/groups/mrivas/ukbb24983/array-combined/snpnet/penalty.v3.rds`
  - New VEP run (version 101, 2020 Aug). We curated its consequence field and grouped into categories, such as PTVs and PAVs. Please read [this page](/17_annotation/20200912_cal_vep_v101).

| Data source | Csq    | n      | Csq_priority | w    |
|-------------|--------|--------|--------------|------|
| Array       | ptv    | 28321  | 1            | 0.5  |
| Array       | pav    | 89161  | 2            | 0.75 |
| HLA         |        | 362    |              | 0.75 |
| Array       | pcv    | 11282  | 3            | 0.9  |
| Array       | intron | 358439 | 4            | 0.9  |
| Array       | utr    | 7928   | 5            | 0.9  |
| Array       | others | 310295 | 6            | 1    |
| CNVs        |        | 275180 |              | 1    |

### version 2

  - `/oak/stanford/groups/mrivas/ukbb24983/array-combined/snpnet/penalty.v2.rds`
  - HLA allelotypes now has .75

### version 1

  - `/oak/stanford/groups/mrivas/ukbb24983/array-combined/snpnet/penalty.rds`
  - PTVs : 0.5, PAVs : 0.75, others : 1.0
  - It turned out that HLA allelotypes have 1.0

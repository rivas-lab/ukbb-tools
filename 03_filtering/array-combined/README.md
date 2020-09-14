# array-combined dataset

- `both_array_list.sh`: the variants in `cal` dataset is measured in two genotyping arrays.


## penalty factor for snpnet

We have penalty factor files for snpnet.

- version 3
  - `/oak/stanford/groups/mrivas/ukbb24983/array-combined/snpnet/penalty.v3.rds`
  - New VEP run.

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

- version 2
  - `/oak/stanford/groups/mrivas/ukbb24983/array-combined/snpnet/penalty.v2.rds`
  - HLA allelotypes now has .75
- version 1
  - `/oak/stanford/groups/mrivas/ukbb24983/array-combined/snpnet/penalty.rds`
  - PTVs : 0.5, PAVs : 0.75, others : 1.0
  - It turned out that HLA allelotypes have 1.0

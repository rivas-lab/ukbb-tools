# Exome 200k data, variant annotation

Yosuke Tanigawa, 2020/10/25-28

Please also check [`exome_oqfe_2020`](/03_filtering/exome_oqfe_2020) documentation as well.

Note: this analysis is not finalized and documentation is still in progress.

## data location

`/oak/stanford/groups/mrivas/ukbb24983/exome/annotation/20201025_exome_oqfe_2020/`

- `UKBexomeOQFE.vep101.tsv.gz`: preliminary version of variant annotation file
- `ukb24983_exomeOQFE.afreq_hwe.20201025.pvar.zst`: allele frequency, missingness, HWE p-value across populations
  - `ukb24983_exomeOQFE.afreq_hwe.20201025.compact.pvar.zst`: subset of fields
- `UKBexomeOQFE.hg19.tsv.gz`: the variants on hg19 coordinates based on liftOver 
- `UKBexomeOQFE.hg19.unmapped.txt.gz`: the list of unmapped variants (based on liftOver)

## Input set of variants

- 17,981,897 variants in the pvar/bim file from [BIM index files for OQFE exome data (Resource 200)](https://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=200).

## variant annotation with VEP

We still have an issue to run Loftee on GRCh38. We apply VEP without Loftee.

- 17,549,650 variants were annotated

## Allele frequency and HWE test

We have 17,777,950 variants in the combined pgen/pvar.zst/psam file with `MAC>=1` filter (see: [`exome_oqfe_2020`](/03_filtering/exome_oqfe_2020)).

We computed the allele frequency and HWE test statistics for the 17,777,950 variants.

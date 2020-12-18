# The exome 200k dataset

Yosuke Tanigawa, 2020/10/26

## file location

`/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020`

- `download`: (symlink to `/scratch`)
- `ukb24983_exomeOQFE.{pgen,pvar.zst,psam}`

## Related contents

- [Variant annotation](/17_annotation/20201025_exome_oqfe_2020)
- [GWAS analysis](/04_gwas/extras/20201026_exome_gwas_parallel)

## methods summary

- [`1_download.sh`](1_download.sh): wrapper script for `gfetch`, submitted with [`1_download.sbatch.sh`](1_download.sbatch.sh)
- [`2_merge.sh`](2_merge.sh): we combine the per-chromosome genotype files into one file.

## QC

### duplicate detection in exome 200k

In the exome 200k dataset in UK Biobank, we find 163 pairs (326 variant IDs) of "duplicated" variants. Each pair represents the exact same variant, if you take `uniq` on (CHROM, POS, REF, ALT) tupple, but they are stored as different variants.

The list of duplicated variants are stored in `/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/ukb24983_exomeOQFE.duplicates.tsv.gz`.

Please see the following two analysis scripts/notebook for more info.

- [`20201217_duplicated_variants_step1.sh`](20201217_duplicated_variants_step1.sh)
- [`20201217_duplicated_variants_step2.ipynb`](20201217_duplicated_variants_step2.ipynb)

### preparing the bim file

```
wget -nd  biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/UKBexomeOQFEbim.zip
```

### `gfetch` software

```
wget -nd  biobank.ctsu.ox.ac.uk/showcase/util/gfetch
chmor 770 gfetch
```



# The exome 200k dataset

Yosuke Tanigawa, 2020/10/26

## file location

`/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020`

- `download`: (symlink to `/scratch`)
- `ukb24983_exomeOQFE.{pgen,pvar.zst,psam}`

## Related contents

- [Variant annotation](/17_annotation/20201025_exome_oqfe_2020)
- [GWAS analysis](/04_gwas/extras/20201026_exome_gwas_parallel)

## QC

The list of QC-passed variants are in `/oak/stanford/groups/mrivas/ukbb24983/exome/qc/oqfe_2020/ukb24983_exomeOQFE.passQC.20201222.lst`

### (cf.) QC criteria in Szustakowski et al 2020

- [Szustakowski et al 2020](https://doi.org/10.1101/2020.11.02.20222232)
  - individual and variant missingness <10%
  - Hardy Weinberg Equilibrium p-value>10^-15
  - minimum read coverage depth of 7 for SNPs and 10 for indels
  - at least one sample per site passed the allele balance threshold > 0.15 for SNPs and 0.20 for indels.

### duplicate detection in exome 200k

In the exome 200k dataset in UK Biobank, we find 163 pairs (326 variant IDs) of "duplicated" variants. Each pair represents the exact same variant, if you take `uniq` on (CHROM, POS, REF, ALT) tupple, but they are stored as different variants.

The list of duplicated variants are stored in `/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/ukb24983_exomeOQFE.duplicates.tsv.gz`.

Please see the following two analysis scripts/notebook for more info.

- [`20201217_duplicated_variants_step1.sh`](20201217_duplicated_variants_step1.sh)
- [`20201217_duplicated_variants_step2.ipynb`](20201217_duplicated_variants_step2.ipynb)

### per-variant missingness and HWE

We performed variant annotation using [VEP](/17_annotation/20201025_exome_oqfe_2020).

### variant-level QC criteria

#### version 2020/12/22

`/oak/stanford/groups/mrivas/ukbb24983/exome/qc/oqfe_2020/ukb24983_exomeOQFE.passQC.20201222.tsv.gz`

The variants passed the following criteria are included: variant-level missingness < 10%, the Hardy-Weinberg equilibrium test (computed within unrelated individuals of white Brisith ancestry) p-value > 10^-15, and the variant is uniqly represented (the CHROM-POS-REF-ALT tupple is uniqly identified) in the PLINK dataset file.

![variant.QC.20201222.UpSetR.png](UpSetR plot summarizing the variant-level QC filter)

**Fig. Summary of variant QC (version 2020/12/22).** In total, we removed 195,920 variants that does not meet any of the following criteria: variant-level missingness < 10% ("missingness"), the Hardy-Weinberg equilibrium test (computed within unrelated individuals of white Brisith ancestry) p-value > 10^-15 ("HWE"), and the variant is uniqly represented (the CHROM-POS-REF-ALT tupple is uniqly identified) in the PLINK dataset file ("duplicated").

## methods summary

- [`1_download.sh`](1_download.sh): wrapper script for `gfetch`, submitted with [`1_download.sbatch.sh`](1_download.sbatch.sh)
- [`2_merge.sh`](2_merge.sh): we combine the per-chromosome genotype files into one file.

### preparing the bim file

```
wget -nd  biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/UKBexomeOQFEbim.zip
```

### `gfetch` software

```
wget -nd  biobank.ctsu.ox.ac.uk/showcase/util/gfetch
chmor 770 gfetch
```

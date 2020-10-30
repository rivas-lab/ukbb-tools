# Exome 200k data, variant annotation

Yosuke Tanigawa, 2020/10/25-30

Please also check [`exome_oqfe_2020`](/03_filtering/exome_oqfe_2020) documentation as well.

## data location

- `/oak/stanford/groups/mrivas/ukbb24983/exome/annotation/20201025_exome_oqfe_2020/`
  - `ukb24983_exomeOQFE.annotation.tsv`: the full variant annotation table.
  - `ukb24983_exomeOQFE.annotation.compact.tsv`: the table with subset of columns.

### column descriptor

The following 20 columns are present in both files.

- `CHROM`: the chromosome
- `POS`: the position
- `ID`: the variant ID
- `REF`: the reference allele
- `ALT`: the alternate allele
- `FILTER`: variant QC filter from VEP - it should be all '.' for now
- `Allele`: the allele used in VEP/Loftee plugin (typically, ALT allele, but it can be different)
- `Csq`: the aggregated consequence group
- `Consequence`: the consequence field from VEP/Loftee
- `SYMBOL`: the gene symbol
- `Gene`: the ensembl ID for the gene
- `f_miss`: the missigness
- `UKB_white_british_hwe_p`: the HWE test p-value
- `UKB_white_british_AF`: the alternate allele frequency computed in white British unrelated individuals in UKB. Note: this is AF, not MAF.
- `UKB_AF`: the alternate allele frequency computed in 200k Exome cohort in UKB. Note: this is AF, not MAF.
- `CHROM`: the chromosome in hg19 (UCSC liftOver)
- `POS`: the position in hg19 (UCSC liftOver)
- `REF`: the reference allele in hg19 (UCSC liftOver)
- `ALT`: the alternate allele in hg19 (UCSC liftOver)
- `liftOver_unmapped_reason`: the error message from UCSC liftOver

The next 133 columns (except `HGVSp`) are present only in the full table.
Among the 133 columns, 65 fields are from VEP/Loftee pipeline.

- `IMPACT`
- `Feature_type`
- `Feature`
- `BIOTYPE`
- `EXON`
- `INTRON`
- `HGVSc`
- `HGVSp`
- `cDNA_position`
- `CDS_position`
- `Protein_position`
- `Amino_acids`
- `Codons`
- `Existing_variation`
- `ALLELE_NUM`
- `DISTANCE`
- `STRAND`
- `FLAGS`
- `VARIANT_CLASS`
- `SYMBOL_SOURCE`
- `HGNC_ID`
- `CANONICAL`
- `MANE`
- `TSL`
- `APPRIS`
- `CCDS`
- `ENSP`
- `SWISSPROT`
- `TREMBL`
- `UNIPARC`
- `GENE_PHENO`
- `SIFT`
- `PolyPhen`
- `DOMAINS`
- `miRNA`
- `HGVS_OFFSET`
- `AF`
- `AFR_AF`
- `AMR_AF`
- `EAS_AF`
- `EUR_AF`
- `SAS_AF`
- `AA_AF`
- `EA_AF`
- `gnomAD_AF`
- `gnomAD_AFR_AF`
- `gnomAD_AMR_AF`
- `gnomAD_ASJ_AF`
- `gnomAD_EAS_AF`
- `gnomAD_FIN_AF`
- `gnomAD_NFE_AF`
- `gnomAD_OTH_AF`
- `gnomAD_SAS_AF`
- `MAX_AF`
- `MAX_AF_POPS`
- `CLIN_SIG`
- `SOMATIC`
- `PHENO`
- `PUBMED`
- `VAR_SYNONYMS`
- `MOTIF_NAME`
- `MOTIF_POS`
- `HIGH_INF_POS`
- `MOTIF_SCORE_CHANGE`
- `TRANSCRIPTION_FACTORS`

Then, we have allele frequency/counts for other UKB populations.

Please see [the annotation file documentation for the array combined dataset](/17_annotation/20201012_array-combined) for column description.

## Input set of variants

- 17,981,897 variants in the pvar/bim file from [BIM index files for OQFE exome data (Resource 200)](https://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=200).
- When preparing the combined plink genotype file, we applied `MAC>=1` filter, which resulted in 17,777,950 variants.


## variant annotation with VEP

We still have an issue to run Loftee on GRCh38. We apply VEP without Loftee.

- 17,549,650 variants were annotated
- `UKBexomeOQFE.vep101.tsv.gz`: preliminary version of variant annotation file

## liftOver (under) to hg19

We apply our liftOver script to find the coordinates on hg19.

- `UKBexomeOQFE.hg19.tsv.gz`: the variants on hg19 coordinates based on liftOver. This file is generated based on the information in the following two files:
  - `UKBexomeOQFE.hg19.mapped.txt.gz`: the list of mapped variants (based on liftOver)
  - `UKBexomeOQFE.hg19.unmapped.txt.gz`: the list of unmapped variants (based on liftOver)

## Allele frequency and HWE test

We have 17,777,950 variants in the combined pgen/pvar.zst/psam file with `MAC>=1` filter (see: [`exome_oqfe_2020`](/03_filtering/exome_oqfe_2020)).

We computed the allele frequency and HWE test statistics for the 17,777,950 variants.

- `ukb24983_exomeOQFE.afreq_hwe.20201025.pvar.zst`: allele frequency, missingness, HWE p-value across populations
  - `ukb24983_exomeOQFE.afreq_hwe.20201025.compact.pvar.zst`: subset of fields


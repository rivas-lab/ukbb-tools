# Variant annotation

## Scripts

- `add_fields.py`: Adds fields for pipe-delimited consequence string and MAF to final tsv.

- `annotate_bims.sh`: The workhorse of the annotation folder. Uses a build of
   VEP to annotate a given BIM file (tab-separated file with CHROM, POS, REF, ALT, MAF).

- `annotate_imp_bims.sh`: Wrapper for submitting jobs to annotate imputation BIMs.

- `combine_cALL.sh`: Concatenates imputation annotations (Autosome only for now).
   Results are written to `/oak/stanford/groups/mrivas/ukbb24983/imp/annotation/annot.tsv.gz`.

- `convert_to_vcf.py`: Converts a given BIM file into a VCF to be annotated by VEP.

- `generate_finngen_bim.sh`: Generates BIM for FinnGen by taking unique positions from a 
   list of summary statistic files.

- `generate_imp_bims.sh`: Generates BIM for FinnGen by taking unique positions from a  
   list of summary statistic files.

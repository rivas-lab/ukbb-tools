# Variant annotation

This folder contains wrappers to the Variant Effect Predictor version 87 that is installed in `/oak/stanford/groups/mrivas/users/guhan/software/ensembl-vep/`.

## Contents

1. [`add_fields.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/17_annotation/add_fields.py): Adds fields for pipe-delimited consequence string and MAF to final tsv.
- Inputs: 
  - `tsv`: Output from tableized, annotated VCF. See `/oak/stanford/groups/mrivas/public_data/loftee/src/tableize_vcf.py` for tableization script.
  - `vcf`: Output from VEP.
  - `maf`: File containing `CHROM`, `POS`, `REF`, `ALT`, and `MAF` as column names.
  - `out`: Output filename.
- Outputs: Tableized VCF with pipe-delimited consequence string and MAF.
- Example usage: `python add_fields.py tsv vep_output.vcf.gz maf.tsv out.tsv`. *NOTE*: This is called within [`annotate_bims.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/17_annotation/annotate_bims.sh); don't call this script on its own unless absolutely necessary.
2. [`annotate_bims.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/17_annotation/annotate_bims.sh): The workhorse of this folder. Uses a build of VEP to annotate a given BIM file.
 - Inputs: 
   - tab-separated file with CHROM, POS, REF, ALT, MAF as the columns. Must have a single-dot extension (e.g. `.txt`, `.tsv`, and not `.txt.tsv`).
   - Genome build (either `GRCh37` or `GRCh38`)
 - Outputs: VEP-annotated `.tsv`.
 - Minimal usage: `sbatch annotate_bims.sh <BIM input filename> <assembly name (GRCh37 or GRCh38)>"`
 - Suggested usage: `sbatch --mem=64000 -t 1-00:00:00 -p mrivas --nodes=1 --cores=8 annotate_bims.sh <BIM input filename> <assembly name (GRCh37 or GRCh38)>"`
 - Example usage: `sbatch --mem=64000 -t 1-00:00:00 -p mrivas --nodes=1 --cores=8 annotate_bims.sh finngen.bim GRCh37`
3. [`annotate_imp_bims.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/17_annotation/annotate_imp_bims.sh): Wrapper for submitting jobs to annotate imputation BIMs.
- Inputs: None.
- Outputs: `/oak/stanford/groups/mrivas/ukbb24983/imp/annotation/imp_annots/ukb24983_imp_chr*_v3.pvar.zst_cf_vep.tsv`. These results are concatenated via [`combine_cALL.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/17_annotation/combine_cALL.sh).
- Example usage: `bash annotate_imp_bims.sh`
4. [`combine_cALL.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/17_annotation/combine_cALL.sh): Concatenates imputation annotation results into a master file, including XY. 
- Inputs: None.
- Outputs: Results are written to `/oak/stanford/groups/mrivas/ukbb24983/imp/annotation/annot.tsv.gz`.
- Example usage: `sbatch combine_cALL.sh`
5. [`convert_to_vcf.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/17_annotation/convert_to_vcf.py): Converts a given BIM file into a VCF to be annotated by VEP.
- Inputs:
  - `.tsv`: File containing `CHROM`, `POS`, `REF`, `ALT`, and `MAF` as column names.
  - `out`: Output filename.
- Outputs: VCF-format file for input into VEP.
- Example usage: `python convert_to_vcf.py $input $temp_vcf`. *NOTE*: This is called within [`annotate_bims.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/17_annotation/annotate_bims.sh); don't call this script on its own unless absolutely necessary.
6. [`generate_finngen_bim.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/17_annotation/generate_finngen_bim.sh): Generates BIM for FinnGen by taking unique positions from a list of summary statistic files.
- Inputs: None.
- Outputs: Unique list of variants (`CHROM`, `POS`, `REF`, `ALT`) in FINNGEN summary statistics (which were once) at `/scratch/groups/mrivas/users/ytanigaw/20200114_FinnGen_R2/summary_stats_hg19/`.
- Example usage: `sbatch generate_finngen_bims.sh`
7. [`generate_imp_bims.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/17_annotation/generate_imp_bims.sh): Generates BIM for imputation data by taking unique positions from a list of summary statistic files.
- Inputs: None.
- Outputs: Unique list of variants (`CHROM`, `POS`, `REF`, `ALT`) in imputation summary statistics at `/oak/stanford/groups/mrivas/ukbb24983/imp/pgen/ukb24983_imp_chr*_v3.pvar.zst`.
- Example usage: `sbatch generate_imp_bims.sh`

## Location of relevant annotation files on `$OAK`

- `cal`: `/oak/stanford/groups/mrivas/private_data/ukbb/variant_filtering/variant_filter_table.tsv.gz`
- `exome`: `/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/spb/data/ukb_exm_spb-[population]-variant_annots_gbe.tsv.gz`
- `imp`: `/oak/stanford/groups/mrivas/ukbb24983/imp/annotation/annot.tsv.gz`

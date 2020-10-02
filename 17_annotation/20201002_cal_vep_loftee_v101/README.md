# The genotyped variant (`cal` dataset) annotation using VEP v101 with loftee plugin
## Yosuke Tanigawa, 2020/10/2

## Data location

- `/oak/stanford/groups/mrivas/ukbb24983/cal/annotation_20201002/ukb24983_cal_cALL_v2_hg19.vep101-loftee.Csq.tsv.gz`

### The column descriptors

The table file contains 79 columns for 805426 variants on the array dataset (`/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_cALL_v2_hg19.pvar`).
The file should have the same ordering of the variants as in the input pvar file.

We have the following information in the 79 (= 5 + 73 + 1) columns.

- The first 5 columns are `#CHROM`, `POS`, `ID`, `REF`, and `ALT`, representing the location, variant ID, and reference and alternate allele of the variant.
- The next 73 columns are the variant annnotation information from VEP + loftee. Please see [`VEP_column_descriptors.md`](VEP_column_descriptors.md) for their description.
- The last column `Csq` represents the grouping of the Consequence field. We have 6 levels. The mapping rule between the `Consequence` field and `Csq` grouping is defined in [`VEP_consequence_group.annotated.tsv`](VEP_consequence_group.annotated.tsv). Note this table is copied from [`../annotation_20200912`](../annotation_20200912).
  - `ptv`: protein-truncating variants
  - `pav`: protein-altering variants
  - `pcv`: protein-coding variants
  - `utr`: UTR-region variants
  - `intron`: intronic-region variants
  - `others`: all other (non-coding) variants

## change log

- 2020/10/2: VEP version 101 with loftee.
  - In September run, we did not use loftee plugin, but now it's incorporated into the pipeline.
  - We now includes all the output fields from VEP + loftee pipeline
- 2020/9/12: VEP version 101 with the latest data.

## Annotation with `VEP` with loftee plugin

We [installed VEP with loftee plugin via Docker/Singularity image](https://github.com/rivas-lab/sherlock-modules/tree/master/vep).

Here, [`vep_example_single_thread.sh`](vep_example_single_thread.sh) illustrates the usage of the software.
We subsequently call the vcf to table conversion script, `/oak/stanford/groups/mrivas/software/loftee/src/tableize_vcf.py`.

### annotation of the array dataset (`cal`)

To minimize the wait time of the computation (wall time), we split the input file into small pieces and apply the pipeline using an array job. We split the input pvar file of 805426 variants into 403 chunks, each of which contains at most 2000 variants. We have the following scripts to facilitate the split, job submission, and merge of the results.

- [`1_split_input_pvar.sh`](1_split_input_pvar.sh): split the input pvar file into pieces.
- [`2_run_vep.sh`](2_run_vep.sh): apply VEP + loftee pipeline. This is basically the same as our example usage script, [`vep_example_single_thread.sh`](vep_example_single_thread.sh).
  - [`2_run_vep.sbatch.submit.sh`](2_run_vep.sbatch.submit.sh): this is just a SLURM array job submission script.
- [`3_combine_vep_results.sh`](3_combine_vep_results.sh): we combine the results (for each chunk) into one table
- [`4_add_Csq_field.ipynb`](4_add_Csq_field.ipynb): we fix the chromosome name (`chr1` vs `1`, etc), add `Csq` field, and sort the table using the input `pvar` file.

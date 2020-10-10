# variant annotation for the CNV dataset

Yosuke Tanigawa, 2020/10/9

We annotated the copy number variants with VEP/Loftee.

Given the variant annotation pipeline is not designed for CNV datasets and we don't know the exact location of the CNVs (because they were called from the array genotype data), one need to be careful when interpreting the annotated variants. Still, the variant annotaion should serve as a resource to locate the genes symbols, for example, around the CNVs of interest.

We were not able to annotate 250/275180 variants - we had hard time annotating these (mostly) long alleles with VEP/loftee.

## Data location

- `/oak/stanford/groups/mrivas/ukbb24983/cnv/annotation_20201003/cnv.vep101-loftee.20201009.Csq.tsv.gz`

This table has all the variants in the pvar file for the CNV dataset (`/oak/stanford/groups/mrivas/ukbb24983/cnv/pgen/cnv.pvar`).
They are NOT sorted in the same order as in the pvar file. Instead, they are sorted by the start and the end position of the CNVs (to enable proper indexing with `tabix`).

- `CHROM`: chromosome
- `POS`: the position as in the pvar file
- `ID`: the variant ID
- `POS_s`: the start position of the CNV (parsed from the ID column)
- `POS_e`: the end position of the CNV (parsed from the ID column) 
- `REF`: the reference column in the pvar file - it is always `N`
- `ALT`: the alternate allele in the pvar file - it is always `+`
- `pvar_order`: the sorting order in the pvar file
- `Allele`: the actual effect allele (parsed from the ID column) - `+` represents insertion and `-` represents deletion
- `Consequence`: the predicted consequence from the variant annotation pipeline.
- `Csq`: the summarized consequence field. Please see [`VEP_consequence_group.md`](/17_annotation/VEP_consequence_group.md) for more information.

Please see [`VEP_column_descriptors.md`](/17_annotation/20201002_cal_vep_loftee_v101/VEP_column_descriptors.md) for the other columns.

## script

- [`1_split_input_pvar.sh`](1_split_input_pvar.sh): we split the pvar file into pieces for efficient computation
- [`2_run_vep.sh`](2_run_vep.sh): a wrapper script to annote variants in a pvar file
- [`2_run_vep.sbatch.submit.sh`](2_run_vep.sbatch.submit.sh): SLURM job submission script
- [`3_combine_vep_results.sh`](3_combine_vep_results.sh): combine the results into one file
- [`4_add_Csq_field.ipynb`](4_add_Csq_field.ipynb): a notebook to add `Csq` field
- [`cnv_pvar_fill_seq.sh`](cnv_pvar_fill_seq.sh): a helper script to fetch actual sequence from fasta file
- [`simplify_cnv_tsv.R`](simplify_cnv_tsv.R): a post processing script to remove redundant information (in REF, ALT, HGVSc, and Allele column)
- [`error_analysis`](error_analysis): some deletion alleles caused errors in variant annotation. we performed error analysis and resubmitted jobs.

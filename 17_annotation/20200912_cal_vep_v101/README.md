# variant annotation on the coding dataset using VEP v101.

## Data location

- `/oak/stanford/groups/mrivas/ukbb24983/cal/annotation_20200912`
- `ukb24983_cal_cALL_v2_hg19.vep101.noLoF_summary.html`: summary fro VEP
- `ukb24983_cal_cALL_v2_hg19.vep101.noLoF.tsv.gz` full output from VEP
- `ukb24983_cal_cALL_v2_hg19.vep101.noLoF.Csq.tsv.gz`: the file with consequence assignment

## Annotation with `VEP`

We used the VEP version 101 (via Docker/Singularity) and annotated the array-variant dataset.

Please also see [`future_development_goals.md`](future_development_goals.md).

## Consequence field

To aggregate `Consequence` field to `Csq` group, we exported the observed set of the predicted `Consequence` and annotated the aggregation rule on Google spreadsheet.

We chcked the definition of each consequence (by looking at the Sequence Ontology [SO]) and grouped them into the following 6 levels.

1. `ptv`: Protein-truncating variant
2. `pav`: Protein-altering variant
3. `pcv`: Protein-coding variant
4. `intron`: Intronic variant
5. `utr`: UTR variant
6. `others`: The others

| Consequence                        | Csq    | n       |
|------------------------------------|--------|---------|
| intron_variant                     | intron | 1646954 |
| downstream_gene_variant            | others | 456023  |
| non_coding_transcript_variant      | others | 416197  |
| upstream_gene_variant              | others | 376652  |
| missense_variant                   | pav    | 266290  |
| intergenic_variant                 | others | 242373  |
| NMD_transcript_variant             | others | 185073  |
| regulatory_region_variant          | others | 171937  |
| non_coding_transcript_exon_variant | others | 116161  |
| 3_prime_UTR_variant                | utr    | 55795   |
| synonymous_variant                 | pcv    | 48269   |
| frameshift_variant                 | ptv    | 47193   |
| splice_region_variant              | pav    | 45968   |
| TF_binding_site_variant            | others | 40215   |
| stop_gained                        | ptv    | 29142   |
| 5_prime_UTR_variant                | utr    | 15676   |
| splice_donor_variant               | ptv    | 10717   |
| splice_acceptor_variant            | ptv    | 8580    |
| inframe_deletion                   | pav    | 1721    |
| start_lost                         | ptv    | 772     |
| inframe_insertion                  | pav    | 729     |
| stop_lost                          | ptv    | 429     |
| coding_sequence_variant            | pcv    | 373     |
| incomplete_terminal_codon_variant  | pcv    | 93      |
| stop_retained_variant              | pcv    | 72      |
| protein_altering_variant           | pav    | 23      |
| start_retained_variant             | pcv    | 14      |
| mature_miRNA_variant               | others | 13      |
| TFBS_ablation                      | others | 12      |

Please see [Google Spreadsheet](https://docs.google.com/spreadsheets/d/11o0Pu7ksyHOS-bhU_ViI8hBfIHWfNDKWH1o3ipHlQ-c/edit?usp=sharing) for more details.

## Reference

- https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html

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


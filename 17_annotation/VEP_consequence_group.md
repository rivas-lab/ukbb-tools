# VEP consequence field grouping

We aggregated `Consequence` field (from `VEP`) to `Csq` groups. We currently have the following 6 levels:

1. `ptv`: Protein-truncating variant
2. `pav`: Protein-altering variant
3. `pcv`: Protein-coding variant
4. `intron`: Intronic variant
5. `utr`: UTR variant
6. `others`: The others

The results of the grouping is written in [`VEP_consequence_group.tsv`](VEP_consequence_group.tsv).

## version history

- 2020/10/9: add `transcript_ablation` and `regulatory_region_ablation` for CNV dataset.
- 2020/9/25: initial grouping: [`20200912_cal_vep_v101` (commit: `922b333d53`)](https://github.com/rivas-lab/ukbb-tools/tree/922b333d53a95731e9d46dd13accf235594f838c/17_annotation/20200912_cal_vep_v101)

## Reference

- https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html
- [The Sequence Ontology (SO)](http://www.sequenceontology.org/)
  - One can download the `obo` file from [the GitHub repo](https://github.com/The-Sequence-Ontology/SO-Ontologies).

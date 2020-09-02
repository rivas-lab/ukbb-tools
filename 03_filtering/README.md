# Filtering

We do not employ sample- and variant-level QC in the traditional sense in this pipeline. Sample-level QC is often done by the BioBank itself before the data is sent out, and variant-level QC is done after GWAS is run (in case any signal can be recovered despite bad data).

The notebooks in this folder do two distinct tasks:

1. [Go over the overall marker quality of different BioBank arrays](https://github.com/rivas-lab/ukbb-tools/blob/master/03_filtering/Marker_QC.ipynb) (marker_QC.ipynb)
2. sample-level QC, define populations, and generate the GWAS covariate file
  - [`population_stratification_20200828`](population_stratification_20200828): has the latest version of our population definition.

## Variant-level QC

We document variant annotation in [17_annotation](/17_annotation) directory.

### array dataset

The variant QC procedure for the genotyped array dataset is well described in the publications from the lab.

1. [DeBoever, C. et al. Medical relevance of protein-truncating variants across 337,205 individuals in the UK Biobank study. Nature Communications 9, 1612 (2018).](https://doi.org/10.1038/s41467-018-03910-9)
2. [Tanigawa, Y. et al. Components of genetic associations across 2,138 phenotypes in the UK Biobank highlight adipocyte biology. Nat Commun 10, 1–14 (2019).](https://doi.org/10.1038/s41467-019-11953-9)

The second paper has the variant QC flowchart as Supplementary Figure 1.

### phased array dataset

Please check [`hap`](hap) directory for more information. As of 3/27/2020, we are asking some issues in BGEN files to UK Biobank.

### imputation dataset

We have [`imp`](imp) directory to document the QC and filtering procedure for the imputation dataset. Specifically, [`imp/2_var_QC`](imp/2_var_QC) has the variant QC procedure.

## Sample-level QC

### version history

- [`population_stratification_20200828`](population_stratification_20200828): update to the most recent remove file; redefine s_asian and e_asian; update local PCs
- [`population_stratification_w24983_20200522.md`](population_stratification_w24983_20200522.md): additional population-specific PCA analysis for White British and "others" cohorts. The "others" cohorts are defined as the set of unrelated individuals who are not captured in the population assignments (based on the Global PC threshold) and defined in [`sample_qc_v3.2.1_pop_count.ipynb`](sample_qc_v3.2.1_pop_count.ipynb). The updated results are also reflected in master phe file, `/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.20200522.phe`. This notebook ([`sample_qc_v3.2.2_pop_covar_update.ipynb`](sample_qc_v3.2.2_pop_covar_update.ipynb)) reads the PCA computation and udpate the sqc file, the GWAS covariate file, and the master phe file.
- `population_stratification_w24983_20200313`: this one corresponds to `sample_qc_v3.2`. The difference between v3.1 and v3.2 are the participant Withdrawal information. We also improved the clarity of the documentation.

### Lists of individuals and variants on Axiom and BiLEVE Arrays

```{bash}
999  awk '($144==1){print $1,$2}' /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200313/ukb24983_master_sqc.20200313.phe > /oak/stanford/groups/mrivas/ukbb24983/sqc/axiom_individuals.txt

1000  awk '($144==0){print $1,$2}' /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200313/ukb24983_master_sqc.20200313.phe > /oak/stanford/groups/mrivas/ukbb24983/sqc/bileve_individuals.txt

 [magu@sh02-ln02 login /oak/stanford/groups/mrivas/private_data/ukbb/variant_filtering]$ zcat variant_filter_table.tsv.gz | awk '($21 ~ /T/){print $5}' > /oak/stanford/groups/mrivas/ukbb24983/sqc/axiom_specific_variants.txt

[magu@sh02-ln02 login /oak/stanford/groups/mrivas/private_data/ukbb/variant_filtering]$ zcat variant_filter_table.tsv.gz | awk '($22 ~ /T/){print $5}' > /oak/stanford/groups/mrivas/ukbb24983/sqc/bileve_specific_variants.txt
```

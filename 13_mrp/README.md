# Multiple Rare-variants and Phenotypes

This subdirectory has all the code necessary to generate the data and figures required for the Multiple Rare-variants and Phenotypes method. The method was first developed into a [bioRxiv manuscript](https://www.biorxiv.org/content/10.1101/257162v5) by DeBoever *et. al.*

## Contents

### Scripts

1. [`mrp_production.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_production.py)
The main script in this folder.
- Required inputs:
  - `--file`: Path to tab-separated file containing list of summary statistic file paths, corresponding studies, phenotypes, and whether or not to use the file in `R_phen` (matrix of correlations across included phenotypes) generation. An example can be found in the [below section](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp#mrp-script-details-and-options).
  - `--metadata_path`: path to tab-separated file containing variants, gene symbols, consequences, MAFs, and LD independence info. An example can be found in the [below section](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp#mrp-script-details-and-options).
  - `--build`: "hg19" or "hg38".
- Outputs: A `log10BF` and corresponding posterior odds for each gene or variant (depending on what is specified in the parameter `M` in  the script) for each parameter specification as described in the [below section](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp#mrp-script-details-and-options).
- Usage/Parameters: See [below section](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp#mrp-script-details-and-options) for details.
- Example usages: See [`mrp_rv_array_biomarkers.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_rv_array_biomarkers.sh), [`mrp_rv_ma_array_bin.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_rv_ma_array_bin.sh), [`mrp_multitrait_array.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_multitrait_array.py), [`mrp_rv_ma_exome_qc_bin.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_rv_ma_exome_qc_bin.sh), and [`mrp_rv_ma_exome_qc_qt.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_rv_ma_exome_qc_qt.sh) for examples of usages and applications.
2. [`mrp_rv_metabolomics.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_rv_metabolomics.sh)
MRP for metabolomics phenotypes without the use of MPC/pLI information.
3. [`mrp_rv_metabolomics_mpc_pli.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_rv_metabolomics_mpc_pli.sh)
MRP for metabolomics phenotypes using MPC/pLI information.
4. [`metabolomics_wrapper.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/metabolomics_wrapper.sh)
Wrapper for [`mrp_rv_metabolomics_mpc_pli.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_rv_metabolomics_mpc_pli.sh).
5. [`combine_mrp_metabolomics.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/combine_mrp_metabolomics.py)
Along with [`combine_mrp_metabolomics.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/combine_mrp_metabolomics.sh), combines all metabolomics results into a single file.
6. [`combine_mrp_metabolomics.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/combine_mrp_metabolomics.sh)
Along with [`combine_mrp_metabolomics.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/combine_mrp_metabolomics.py), combines all metabolomics results into a single file.
7. [`mrp_rv_array_biomarkers.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_rv_array_biomarkers.sh)
MRP for array data for the 35 biomarkers focused on in [Sinnott-Armstrong, et. al.](https://www.nature.com/articles/s41588-020-00757-z).
8. [`mrp_rv_exome_biomarkers.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_rv_exome_biomarkers.sh)
MRP for exome data for the 35 biomarkers focused on in [Sinnott-Armstrong, et. al.](https://www.nature.com/articles/s41588-020-00757-z) using MPC/pLI information.
9. [`mrp_rv_exome_var_biomarkers.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_rv_exome_var_biomarkers.sh)
MRP for exome data for the 35 biomarkers focused on in [Sinnott-Armstrong, et. al.](https://www.nature.com/articles/s41588-020-00757-z) without using MPC/pLI information.
10. [`mrp_rv_ma_array_biomarkers.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_rv_ma_array_biomarkers.sh)
MRP for meta-analysis on array data for the 35 biomarkers focused on in [Sinnott-Armstrong, et. al.](https://www.nature.com/articles/s41588-020-00757-z).
11. [`mrp_rv_ma_array_bin.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_rv_ma_array_bin.sh)
MRP for meta-analysis on array data for all binary traits in the GBE.
12. [`mrp_rv_ma_array_qt.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_rv_ma_array_qt.sh)
MRP for meta-analysis on array data for all quantitative traits in the GBE.
13. [`mrp_rv_ma_exome_qc_bin.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_rv_ma_exome_qc_bin.sh)
MRP for meta-analysis on exome data for all binary traits in the GBE.
14. [`mrp_rv_ma_exome_qc_qt.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_rv_ma_exome_qc_qt.sh)
MRP for meta-analysis on exome data for all quantitative traits in the GBE.
15. [`combine_mrp.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/combine_mrp.py)
Along with [`combine_mrp.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/combine_mrp.sh), combines all MRP results into a single file.
16. [`combine_mrp.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/combine_mrp.sh)
Along with [`combine_mrp.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/combine_mrp.py), combines all MRP results into a single file.
17. [`generate_gbe_tables.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/generate_gbe_tables.sh)
Generates GBE-friendly tables from the results gathered in [`combine_mrp.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/combine_mrp.sh).
18. [`biomarker_ldsc_wrapper.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/biomarker_ldsc_wrapper.sh)
Performs pairwise LD score regressions across the 35 biomarker traits in [Sinnott-Armstrong, et. al.](https://www.nature.com/articles/s41588-020-00757-z).
19. [`biomarkers_rg_agg.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/biomarkers_rg_agg.sh)
Aggregates the results from [`biomarker_ldsc_wrapper.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/biomarker_ldsc_wrapper.sh).
20. [`biomarkers_cluster.R`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/biomarkers_cluster.R)
Generates a distance matrix and clusters the 35 biomarkers from [Sinnott-Armstrong, et. al.](https://www.nature.com/articles/s41588-020-00757-z).
21. [`mrp_multitrait_array.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_multitrait_array.py)
Performs multitrait analysis, one per cluster, on array data.
22. [`mrp_multitrait_exome.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_multitrait_exome.py)
Performs multitrait analysis, one per cluster, on exome data.
23. [`mutplot_input.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mutplot_input.py)
Generates amino-acid information from HGVSp for input into [`alpl.ipynb`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/manhattan.ipynb).

### Notebooks

1. [`alpl.ipynb`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/alpl.ipynb)
Generates Supplementary Figure S4 from [`alpl.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/alpl.tsv).
2. [`array-exome.ipynb`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/array-exome.ipynb)
Generates Supplementary Table S1 ([`exome_array_diff.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/exome_array_diff.tsv)), Supplementary Figure S2, and Supplementary Table S3 ([`mpc_pli_diff.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mpc_pli_diff.tsv)).
3. [`manhattan.ipynb`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/manhattan.ipynb)
Generates the six sub-figures of Figure 3 (see [Images](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp#images) subsection below).
4. [`power_comparison.ipynb`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/power_comparison.ipynb)
Generates Figures 2B and 4C.


### Images

1. [`Bone and Joint.png`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/Bone%20and%20Joint.png)
Manhattan plot from Figure 3 [notebook](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/manhattan.ipynb) for Alkaline phosphatase, Calcium, and Vitamin D.
2. [`Cardiovascular.png`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/Cardiovascular.png)
Manhattan plot from Figure 3 [notebook](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/manhattan.ipynb) for Apolipoprotein A, Apolipoprotein B, C-reactive protein, Cholesterol, HDL cholesterol, LDL cholesterol, Lipoprotein A, and Trigylcerides.
3. [`Diabetes.png`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/Diabetes.png)
Manhattan plot from Figure 3 [notebook](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/manhattan.ipynb) for Glucose and HbA1c.
4. [`Hormone.png`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/Hormone.png)
Manhattan plot from Figure 3 [notebook](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/manhattan.ipynb) for IGF-1, SHBG, and Testosterone.
4. [`Liver.png`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/Liver.png)
Manhattan plot from Figure 3 [notebook](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/manhattan.ipynb) for AST to ALT ratio, Alanine aminotransferase, Albumin, Aspartate aminotransferase, Direct bilirubin, Gamma glutamyltransferase, and Total bilirubin.
5. [`Renal.png`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/Renal.png)
Manhattan plot from Figure 3 [notebook](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/manhattan.ipynb) for Creatinine, Cystatin C, Microalbumin in urine, Non-albumin protein, Phosphate, Total protein, Urate, Urea, and eGFR.

### Other Files

1. [`alpl.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/alpl.tsv)
Captures amino acid information within the ALPL gene (see [`mutplot_input.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mutplot_input.py)). Used as input to make Supplementary Figure S4 (see [`alpl.ipynb`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/alpl.ipynb)).
2. [`biomarkers_clusters.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/biomarkers_clusters.tsv)
Captures biomarker cluster information (see [`biomarkers_cluster.R`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/biomarkers_cluster.R)).
3. [`exome_array_diff.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/exome_array_diff.tsv)
Captures those gene-trait associations in which the jump from array to exome results in a substantial increase in power (Supplementary Table S1 - see [`array-exome.ipynb`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/array-exome.ipynb)).
4. [`gbe_blacklist.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/gbe_blacklist.tsv)
A list of "blacklisted" GBE phenotypes (to exclude from final results). Curated manually.
5. [`metabolomics_phenos.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/metabolomics_phenos.tsv)
A list of metabolomics phenotypes given from Nightingale.
6. [`mpc_pli_diff.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mpc_pli_diff.tsv)
Captures those gene-trait associations in which the incorporation of MPC and pLI results in a substantial increase in power (Supplementary Table S3 - see [`array-exome.ipynb`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/array-exome.ipynb)).
7. [`pav_nasa.txt`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/pav_nasa.txt)
A list of protein-altering variant assocations found in [Sinnott-Armstrong, et. al.](https://www.nature.com/articles/s41588-020-00757-z).
8. [`ptv_nasa.txt`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/ptv_nasa.txt)
A list of protein-truncating variant assocations found in [Sinnott-Armstrong, et. al.](https://www.nature.com/articles/s41588-020-00757-z).
9. [`requirements.txt`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/requirements.txt)
A list of package requirements for MRP.
10. [`sumstat_paths.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/sumstat_paths.tsv)
Contains summary statistic paths, among other metadata, for the 35 biomarkers.

## Pipelines and Workflows

## MRP: Script details and options

A full list of options can be obtained by running `python3 mrp_production.py -h` from the command line. This output is replicated here for reference:

```{bash}
(base) user$ python mrp_production.py -h
usage: mrp_production.py [-h] --file MAP_FILE --metadata_path METADATA_PATH
                         --build {hg19,hg38}
                         [--chrom {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y} [{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y} ...]]
                         [--mean MEAN]
                         [--R_study {independent,similar} [{independent,similar} ...]]
                         [--R_var {independent,similar} [{independent,similar} ...]]
                         [--M {variant,gene} [{variant,gene} ...]]
                         [--sigma_m_types {sigma_m_mpc_pli,sigma_m_var,sigma_m_1,sigma_m_005} [{sigma_m_mpc_pli,sigma_m_var,sigma_m_1,sigma_m_005} ...]]
                         [--variants {pcv,pav,ptv,all} [{pcv,pav,ptv,all} ...]]
                         [--maf_thresh MAF_THRESHES [MAF_THRESHES ...]]
                         [--se_thresh SE_THRESHES [SE_THRESHES ...]]
                         [--prior_odds PRIOR_ODDS_LIST [PRIOR_ODDS_LIST ...]]
                         [--p_value {farebrother,davies,imhof} [{farebrother,davies,imhof} ...]]
                         [--exclude EXCLUDE] [--filter_ld_indep]
                         [--out_folder OUT_FOLDER]
                         [--out_filename OUT_FILENAME]

MRP takes in several variables that affect how it runs.

optional arguments:
  -h, --help            show this help message and exit
  --file MAP_FILE       path to tab-separated file containing list of: 
                                 summary statistic file paths,
                                 corresponding studies,
                                 phenotypes, and
                                 whether or not to use the file in R_phen generation.
                               
                                 format:
                                 
                                 path        study        pheno        R_phen
                                 /path/to/file1   study1    pheno1     TRUE
                                 /path/to/file2   study2    pheno1     FALSE
                                 
  --metadata_path METADATA_PATH
                        path to tab-separated file containing:
                                 variants,
                                 gene symbols,
                                 consequences,
                                 MAFs,
                                 and LD independence info.
                               
                                 format:
                                 
                                 V       gene_symbol     most_severe_consequence maf  ld_indep
                                 1:69081:G:C     OR4F5   5_prime_UTR_variant     0.000189471     False
                                
  --build {hg19,hg38}   genome build (hg19 or hg38). Required.
  --chrom {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y} [{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y} ...]
                        chromosome filter. options include 1-22, X, and Y
  --mean MEAN           prior mean of genetic effects (Default: 0).
  --R_study {independent,similar} [{independent,similar} ...]
                        type of model across studies. 
                                 options: independent, similar (default: similar). can run both.
  --R_var {independent,similar} [{independent,similar} ...]
                        type(s) of model across variants. 
                                 options: independent, similar (default: independent). can run both.
  --M {variant,gene} [{variant,gene} ...]
                        unit(s) of aggregation. 
                                 options: variant, gene (default: gene). can run both.
  --sigma_m_types {sigma_m_mpc_pli,sigma_m_var,sigma_m_1,sigma_m_005} [{sigma_m_mpc_pli,sigma_m_var,sigma_m_1,sigma_m_005} ...]
                        scaling factor(s) for variants.
                                 options: var (i.e. 0.2 for ptvs, 0.05 for pavs/pcvs), 
                                 1, 0.05 (default: mpc_pli). can run multiple.
  --variants {pcv,pav,ptv,all} [{pcv,pav,ptv,all} ...]
                        variant set(s) to consider. 
                                 options: proximal coding [pcv], 
                                          protein-altering [pav], 
                                          protein truncating [ptv],
                                          all variants [all]
                                          (default: ptv). can run multiple.
  --maf_thresh MAF_THRESHES [MAF_THRESHES ...]
                        which MAF threshold(s) to use. must be valid floats between 0 and 1 
                                 (default: 0.01).
  --se_thresh SE_THRESHES [SE_THRESHES ...]
                        which SE threshold(s) to use. must be valid floats between 0 and 1 
                                 (default: 0.2). NOTE: This strict default threshold is best suited for binary
                                 summary statistics. For quantitative traits, we suggest the use of a higher
                                 threshold.
  --prior_odds PRIOR_ODDS_LIST [PRIOR_ODDS_LIST ...]
                        which prior odds (can be multiple) to use in calculating posterior 
                                 probabilities. must be valid floats between 0 and 1 (default: 0.0005, expect 
                                 1 in 2000 genes to be a discovery).
  --p_value {farebrother,davies,imhof} [{farebrother,davies,imhof} ...]
                        which method(s) to use to convert Bayes Factors to p-values. if command 
                                 line argument is invoked but method is not specified, will throw an error 
                                 (i.e., specify a method when it is invoked). if not invoked, p-values will not 
                                 be calculated. options: farebrother, davies, imhof. NOTE: --p_value imports R 
                                 objects and methods, which slows down MRP. farebrother is fastest and 
                                 recommended if p-values are a must.
  --exclude EXCLUDE     path to file containing list of variants to exclude from analysis.
                        
                                 format of file:
                        
                                 1:69081:G:C
                                 1:70001:G:A
                                
  --filter_ld_indep     whether or not only ld-independent variants should be kept (default: False;
                                 i.e., use everything).
  --out_folder OUT_FOLDER
                        folder to which output(s) will be written (default: current folder).
                                 if folder does not exist, it will be created.
  --out_filename OUT_FILENAME
                        file prefix with which output(s) will be written (default: underscore-delimited
                                 phenotypes).
```

## MRP: Variant groupings

The following groups are how we assign priors (`sigma_m`) to variants.

`ptv = ['frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'start_lost', 'stop_lost']`

`pav = ['protein_altering_variant', 'inframe_deletion', 'inframe_insertion', 'splice_region_variant', 'start_retained_variant', 'stop_retained_variant', 'missense_variant']`

`proximal_coding = ['synonymous_variant', '5_prime_UTR_variant', '3_prime_UTR_variant', 'coding_sequence_variant', 'incomplete_terminal_codon_variant', 'TF_binding_site_variant']`

`intron = ['regulatory_region_variant', 'intron_variant', 'intergenic_variant', 'downstream_gene_variant', 'mature_miRNA_variant', 'non_coding_transcript_exon_variant', 'upstream_gene_variant', 'NA', 'NMD_transcript_variant']`

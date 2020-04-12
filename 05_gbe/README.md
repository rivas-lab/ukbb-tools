# Processing data for the Global Biobank Engine — a pipeline

## Contents

1. [`combine_phe.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/combine_phe.py): Makes a new version of a `master.phe` file based on the contents of [`phenotype_info.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/phenotype_info.tsv).
- Inputs: None.
- Outputs: In order to run PheWAS, we need a master table containing all the phenotype data we've processed (see the [Phenotyping](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/) folder for more info). This maker script will walk through the [`phenotype_info.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/phenotype_info.tsv) file and select phenotypes with `N_GBE` (number of white British individuals) >= 100 for inclusion. It then loads each of those files into memory (!!) via a dictionary keyed on sample ID, then writes all the data out to file. A new `master.YYYYMMDD.phe` and corresponding `.info.tsv` file (the latter of which is a strict subset of the [`phenotype_info.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/phenotype_info.tsv) file) will be created at `/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/`. 
- Example usage: `sbatch -p normal,mrivas,owners --mem=128000 -t 1-00:00:00 -J masterPhe --wrap="python ../../05_gbe/combine_phe.py"`. *NOTE*: [`combine_phe.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/combine_phe.py) is called by the [`02_phenotyping/scripts/update_phe_icd_master.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/update_phe_icd_master.sh) wrapper script, if one is interested in running all of the GBE-info-related scripts together.
2. [`exome_phenotype_info.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/exome_phenotype_info.tsv):  Contains fields, tables, baskets, N overall and per major population (white British, non-British white, African, East Asian, South Asian), source `.tsv` files, dates that the info files were generated, and paths to the phenotypes for exome data.
3. `icdinfo.[population].shortnames.txt`: Contain GBE ID, N (of the population in the filename), GBE names, and GBE "clean" names for each defined phenotype.
4. [`icdinfo.shortnames.BBJ.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/icdinfo.shortnames.BBJ.tsv): Contains GBE ID, N (in Biobank Japan), GBE names, and GBE "clean names" for each phenotype in Biobank Japan that we have mapped back to our own.
5. [`icdinfo.shortnames.MVP.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/icdinfo.shortnames.MVP.tsv): Contains GBE ID, N (in the Million Veterans Program), GBE names, and GBE "clean names" for each phenotype in the Million Veterans Program that we have mapped back to our own.
6. [`icdinfo.shortnames.exome.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/icdinfo.shortnames.exome.tsv): Contains GBE ID, N (in the white British population for exome-sequenced individuals), GBE names, and GBE "clean names" for each defined phenotype.
7. [`icdinfo.shortnames.generator.ipynb`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/icdinfo.shortnames.generator.ipynb): Generates the initial `icdinfo.shortnames.tsv` that has since been moved to [this Google Sheet of "clean" GBE names](http://bit.ly/GBE_names).
8. [`icdinfo.shortnames.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/icdinfo.shortnames.tsv): Contains GBE ID, N (in the white British population), GBE names, and GBE "clean names" for each defined phenotype.
9. [`icdinfo.txt`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/icdinfo.white_british.txt): A symlink to [`icdinfo.white_british.txt`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/icdinfo.white_british.txt), which is a GBE index file. This is a strict subset of the information contained in [`phenotype_info.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/phenotype_info.tsv).
10. [`icdinfo_with_shortnames.ipynb`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/icdinfo_with_shortnames.ipynb): Generates population-specific `icdinfo` files with shortnames as pulled from [this Google Sheet of "clean" GBE names](http://bit.ly/GBE_names).
11. [`make_exome_phenotype_info.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/make_exome_phenotype_info.py): Updates [`exome_phenotype_info.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/exome_phenotype_info.tsv) by intersecting phenotype files with population definitions in `/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/current/info`.
- Inputs: None.
- Outputs: A new [`exome_phenotype_info.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/exome_phenotype_info.tsv) file.
- Example usage: `sbatch -p normal,mrivas,owners --mem=64000 -t 1-00:00:00 -J exm_phen_info --wrap="python ../../05_gbe/make_exome_phenotype_info.py"`. *NOTE*: This is called by the [`02_phenotyping/scripts/update_phe_icd_master.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/update_phe_icd_master.sh) wrapper script, if one is interested in running all of these scripts together.
12. [`phenotype_info.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/phenotype_info.tsv): Contains fields, tables, baskets, N overall and per major population (white British, non-British white, African, East Asian, South Asian), source `.tsv` files, dates that the info files were generated, and paths to the phenotypes for array data.
13. [`update_icdinfo.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/update_icdinfo.sh): Updates [`icdinfo.txt`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/icdinfo.white_british.txt) for a given population (white British by default).
- Inputs: None.
- Outputs: A new [`icdinfo.white_british.txt`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/icdinfo.white_british.txt) file, and symlink at [`icdinfo.txt`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/icdinfo.white_british.txt).
- Example usage: `bash update_icdinfo.sh [population]`. *NOTE*: This is called by the [`02_phenotyping/scripts/update_phe_icd_master.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/update_phe_icd_master.sh) wrapper script, if one is interested in running all of these scripts together.
14. [`update_phe_info.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/update_phe_info.sh): Updates [`phenotype_info.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/phenotype_info.tsv) by concatenating all info files in the symlink directory, `/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/current/info`.
- Inputs: None.
- Outputs: A new [`phenotype_info.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/phenotype_info.tsv) file.
- Example usage: `bash update_phe_info.sh`. *NOTE*: This is called by the [`02_phenotyping/scripts/update_phe_icd_master.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/update_phe_icd_master.sh) wrapper script, if one is interested in running all of these scripts together.

## Pipelines and Workflows

### Generating requisite information files for updating summary statistics

*NOTE:* The [`02_phenotyping/scripts/update_phe_icd_master.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/update_phe_icd_master.sh) wrapper script calls the smaller scripts in this directory, and will update [`phenotype_info.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/phenotype_info.tsv), [`exome_phenotype_info.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/exome_phenotype_info.tsv), [`icdinfo.txt`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/icdinfo.txt) in this directory. It will also create a new `master.YYYYMMDD.phe` and corresponding `.info.tsv` file (a strict subset of the [`phenotype_info.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/phenotype_info.tsv) file) at `/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/` via a batch job. See the [complete pipeline](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping#generating-and-updating-phenotypes-and-summary-statistics) involving some of these scripts. Generally, after this wrapper script is run, you are ready for [GWAS](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas).

# "Clean" GBE names

We maintain a [Google spreadsheet](http://bit.ly/GBE_names) that contains "clean" GBE names (publication-ready for figures, etc.). This is (generally) a more updated version of [`icdinfo.shortnames.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/icdinfo.shortnames.tsv).

# Updating the `icdinfo.[population].shortnames.txt` file across populations

In order to generate population-specific `icdinfo.[population].txt` files, we run the following command:

```{bash}
for pop in white_british non_british_white s_asian e_asian african; do bash update_icdinfo.sh ${pop}; done
```

Then, we run the R notebook, [`icdinfo_with_shortnames.ipynb`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/icdinfo_with_shortnames.ipynb), to add the ["clean" GBE names](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe#clean-GBE-names) to these files, and delete the intermediate files with the following command:

```{bash}
for pop in non_british_white s_asian e_asian african; do rm icdinfo.${pop}.txt; done
```

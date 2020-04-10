# Phenotyping scripts

This folder's scripts help define phenotypes within our pipeline.

## Contents

1. [`annotate_phe.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/annotate_phe.py): Contains an important helper function called `make_phe_info`, which makes `.info` files based on information provided.
- Inputs: None; this helper function is called within [`tsv_to_phenos.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/tsv_to_phenos.py) in order to make `.info` files alongside `.phe` files. In the `main()` of [`annotate_phe.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/annotate_phe.py), the info files are made for all phenotypes in `phenotype_info.tsv`; *this is not recommended*. The main function of this script is to supply the helper function.
- Outputs: `.info` files that are made within the `/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/` directory.
- Example usage: *Not recommended*. See instead [`tsv_to_phenos.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/tsv_to_phenos.py).
2. [`check_for_blank_fields.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/check_for_blank_fields.py): This is meant as a quick sanity check that our annotation table is as complete as can be before a phenotyping session.
- Inputs: None.
- Outputs: None. Prints an informative message regarding data fields that we possess that do not have metadata in the Data Dictionary Showcase.
- Example usage: `python check_for_blank_fields.py`
3. [`make_phe.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/make_phe.py): Contains helper functions `create_bin_phe_file` and `create_qt_phe_file`, which do the grunt work of making phenotypes based on input `.tab` files.
- Inputs: None; these helper functions are called within [`tsv_to_phenos.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/tsv_to_phenos.py) in order to make `.phe` files. In the `main()` of [`make_phe.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/make_phe.py), there is an argument parser, *but using this is not recommended*. The main function of this script is to supply the helper functions.
- Outputs: `.phe` files that are made within the `/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/` directory.
- Example usage: *Not recommended*. See instead [`tsv_to_phenos.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/tsv_to_phenos.py).
4. [`tsv_to_phenos.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/tsv_to_phenos.py): Turns an annotation `.tsv` into lab phenotypes. 
- Inputs: The master annotation table [`tables/ukb_annotations.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/tables/ukb_annotations.tsv) is used by default under-the-hood, but any `.tsv` (that has the same column headers in the same order, ideally) can be specified using `--tsv`. Also, optionally, the flag `--missing-is-control` is supplied if one wants to make phenotypes using missing data as controls for a given annotation.
- Outputs: Phenotype (3-column FID, IID, Phenotype) files to their respective subdirectory within `/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/`. By default, only phenotypes that have not been generated yet will be.
- Example usage: `sbatch -p mrivas -t 01-00:00:00 --mem=64000 --wrap="python tsv_to_phenos.py"`
5. [`update_annotations.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/update_annotations.py): Creates a new [`tables/ukb_annotations.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/tables/ukb_annotations.tsv) file using previous annotations and a Data Dictionary Showcase file.
- Inputs: None.
- Outputs: A new annotation file within `../tables/annotations`, and a symlink, [`tables/ukb_annotations.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/tables/ukb_annotations.tsv). *NOTE*: If a filename already exists with the current date in the `../tables/annotations` folder (e.g. `../tables/annotations/ukb_YYYYMMDD.tsv`), this will be overwritten. In the unlikely event that data refreshes happen multiple times in one day, feel free to rename an older version so there is no namespace conflict.
- Example usage: `python update_annotations.py`
6. [`update_field_table_map.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/update_field_table_map.sh): Creates a new [`tables/field_table_basket_date.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/tables/field_table_basket_date.tsv) file.
- Inputs: None. 
- Outputs: [`tables/field_table_basket_date.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/tables/field_table_basket_date.tsv), which is a map of fields to tables, baskets, and release dates. A field might have multiple baskets/tables; phenotypes that are made using [`tsv_to_phenos.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/tsv_to_phenos.py) take the latest version of a field across baskets and tables by sorting by release date.
- Example usage: `bash update_field_table_map.sh`
7. [`update_phe_icd_master.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/update_phe_icd_master.sh): Wrapper that updates info files and master.phe.
- Inputs: None. 
- Outputs: [`phenotype_info.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/phenotype_info.tsv), [`exome_phenotype_info.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/exome_phenotype_info.tsv), and [`icdinfo.txt`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/icdinfo.txt) reference files in [`05_gbe`](https://github.com/rivas-lab/ukbb-tools/tree/master/05_gbe), and a new `master.phe` and corresponding `.info` file at `/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/`.
- Example usage: `bash update_phe_icd_master.sh`. See the calls to smaller scripts in [`05_gbe`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/).
s

## Pipelines and Workflows

See [parent folder](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping) for the [complete pipeline](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping#generating-and-updating-phenotypes-and-summary-statistics) involving most of these scripts.

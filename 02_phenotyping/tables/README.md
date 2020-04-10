# Phenotyping tables

This is a holding tank for tables relevant to phenotyping.

## Contents

1. [`annotations`](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/tables/annotations): Contains all versions of [`ukb_annotations.tsv`](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/tables/ukb_annotations.tsv).
2. [`Data_Dictionary_Showcases`](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/tables/Data_Dictionary_Showcases): Contains all versions of [`Data_Dictionary_Showcase.csv`](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/tables/Data_Dictionary_Showcase.csv).
3. [`old`](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/tables/old): Old annotation tables that don't have standardized column names or orders.
4. [`Data_Dictionary_Showcase.csv`](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/tables/Data_Dictionary_Showcase.csv): Symlink to the most recent table in [`Data_Dictionary_Showcases`](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/tables/Data_Dictionary_Showcases); the UK Biobank's annotation file.
5. [`field_table_basket_date.tsv`](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/tables/field_table_basket_date.tsv`): A map of fields to tables, baskets, and release dates. Note that a field might have multiple baskets/tables.
6. [`ukb_annotations.tsv`](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/tables/ukb_annotations.tsv): Symlink to the most recent table in [`annotations`](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/tables/annotations); our master annotation file.
7. [`ukbb24983_tables.tsv`](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/tables/ukbb24983_tables.tsv): A list of `.tab` files, their release dates, and their respective baskets.

## Pipelines and Workflows

See [parent folder](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping) for the [complete pipeline](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping#generating-and-updating-phenotypes-and-summary-statistics) involving most of these tables.
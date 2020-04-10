# Phenotyping

The contents of this folder will take you through how to define phenotypes within our pipeline.

## Contents
1. [`extras`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/extras): Scripts for generation and annotation of non-pipeline/"extra" phenotypes. Check README.md within this folder for extras-specific details. 
2. [`scripts`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts): Scripts for generation and annotation of pipeline phenotypes. Check README.md within this folder for script-specific details. 
3. [`tables`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/tables): Relevant reference and generated tables from the pipeline. Check README.md  within this folder for table-specific details. 

## Pipelines and Workflows

### Generating and updating phenotypes and summary statistics

Let's say a new table (or tables) has/have come in from the UK Biobank, and you are ready to turn it into phenotypes for the lab. Follow the following steps:

1. Download and convert the new table(s) using the instructions in the [relevant section in `01_preprocessing`](https://github.com/rivas-lab/ukbb-tools/blob/master/01_preprocessing#unpackingdecryptingconverting-the-data). Make sure that the resultant `.tab` file(s) is/are in the proper place.
2. [Download](http://biobank.ctsu.ox.ac.uk/~bbdatan/Data_Dictionary_Showcase.csv) a new copy of Data_Dictionary_Showcase.csv, rename it as `Data_Dictionary_Showcase.YYYYMMDD.csv`, and place it in [`tables/Data_Dictionary_Showcases`](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/tables/Data_Dictionary_Showcases). Symlink the existing [Data_Dictionary_Showcase.csv`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/tables/Data_Dictionary_Showcase.csv) to it (i.e., run `ln -sf tables/Data_Dictionary_Showcases/Data_Dictionary_Showcase.YYYYMMDD.csv tables/Data_Dictionary_Showcase.csv`).
3. [This Google Drive](http://bit.ly/UKB24983_tables) spreadsheet contains a list of tables and their respective baskets, descriptions, release dates, and paths on Sherlock. Make sure that your new table(s) has/have a row dedicated to it in this spreadsheet. Then, download a new copy of this table, rename it as `ukbb24983_tables.tsv`, and replace the one at [`tables/ukbb24983_tables.tsv`](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/tables/ukbb24983_tables.tsv).
4. Run [`scripts/update_field_table_map.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/update_field_table_map.sh) to generate a new [`tables/field_table_basket_date.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/tables/field_table_basket_date.tsv). As per the name, the resultant file contains a map between UK Biobank fields, their respective tables, and the release dates of those tables.
5. Run [`scripts/check_for_blank_fields.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/check_for_blank_fields.py) to generate a list of fields that are not in [`tables/Data_Dictionary_Showcase.csv`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/tables/Data_Dictionary_Showcase.csv).
- *NOTE*: As of 4/9/2020, 99 of the fields that we have in our data are now considered "defunct" and are thus no longer included in Data Dictionary Showcases. We have chosen to ignore these 99. These 99 will pop up as a result of this script. If you get a number of missing phenotypes that is any more than 99, email access@ukbiobank.ac.uk to get these into a new Data Dictionary Showcase (or otherwise determine whether they are also defunct), then run step 2 again if necessary.
6. Run [`scripts/update_annotations.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/update_annotations.py) to generate a new `tables/ukb_annotations.tsv` file.
7. *Important: If we choose to have a phenotyping session for the new data, this is where it must be done in the pipeline.* Upload a copy of the [`tables/ukb_annotations.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/tables/ukb_annotations.tsv) file to Google Drive, and annotate as described in [the section below](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping#phenotyping-sessions-what-is-a-phenotyping-session-why-do-we-do-it).
8. Run [`scripts/tsv_to_phenos.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/tsv_to_phenos.py) (which, by default, operates on [`tables/ukb_annotations.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/tables/ukb_annotations.tsv) to generate all new/updated `.phe` files.
- *NOTE*: Depending on how many new or updated phenotypes are in the new data table(s), this may take quite some time. For context, generating all ~2300 phenotypes took around 6 days with a job script that requested 64GB of memory on Sherlock.
9. Run [`scripts/update_phe_icd_master.sh`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/scripts/update_phe_icd_master.sh), a wrapper that updates the [`phenotype_info.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/phenotype_info.tsv), [`exome_phenotype_info.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/exome_phenotype_info.tsv), and [`icdinfo.txt`](https://github.com/rivas-lab/ukbb-tools/blob/master/05_gbe/icdinfo.txt) reference files in [`05_gbe`](https://github.com/rivas-lab/ukbb-tools/tree/master/05_gbe) and generates a new `master.phe` and corresponding `.info` file at `/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/`. The separate calls to the different scripts are also easy to make out within this script, if you need to run one or more parts of Step 8 separately. 
10. That's it! You're now ready to [run new GWAS](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas).

### Phenotyping sessions: What is a phenotyping session? Why do we do it?

A phenotyping session is the first step in this process, and should be convened when new data is added to our existing dataset. Decisions made at a phenotyping session allow the lab to extract binary (BIN) or quantitative (INI/QT) phenotypes from new or updated Biobank fields, and subsequently allow the lab to run a number of downstream analyses including GWAS and PheWAS.

#### Getting the group together

This might just be the hardest part of this section - you'll need to find a time when all of the group is ready, willing, and able to work on this together. Doing it by your lonesome is sucky and should never happen. Ping the Slack well in advance and get a quorum of people together before these things happen.

#### At the session

At the beginning of the session, you will need to transfer the most current copy of [`tables/ukb_annotations.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/tables/ukb_annotations.tsv) to a Google Sheet after running steps 1-6 as described [above](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping#generating-and-updating-phenotypes-and-summary-statistics).

At the Biobank phenotyping session, the concept is to divide and conquer. You'll be given a number of phenotypes that don't yet have annotations to define. Look up the field ID [in the Biobank search box](http://biobank.ctsu.ox.ac.uk/crystal/search.cgi). This should give you the name of the field ID when you click on the corresponding number. Fill this name in the appropriate column in the Google sheet. Then, for each row, you should ask yourself:

1) Is this field referring to one trait or multiple? 
- One trait would be something like "Have you smoked before in your life?"
- Multiple traits would be something like "Please mark if you have had stroke, cardiac event, both, or none" (you would be splitting this into two phenotypes - "Had stroke" and "Had cardiac event").
    - In the case of multiple-trait fields, you will need to copy paste that row that you're on as many times as needed (in the case above, copy and paste it once to have a total of 2 rows). Assign indices to each of the rows in the `QT_index` or `BIN_index` columns, and put the total number of traits in the corresponding `QT_total_num` and `BIN_total_num` columns.
  
2) Assuming that the row now refers to a single trait and not multiple, is this a binary trait (BIN) or a quantitative (INI/QT) trait? 
- Depending on whether or not the value is binary or quantitative, assign `BIN_total_num` or `QT_total_num` to be 1 and the other to be 0. This is important for automatically generating GBE IDs, discussed below.
- In the case of a binary trait, your decisions are pretty simple. Find out what value corresponds to cases and what value corresponds to controls via the [Biobank search box](http://biobank.ctsu.ox.ac.uk/crystal/search.cgi).
  - Specifically, search for the field ID of interest and click on it. You should see some tabs, with "Data" being highlighted by default. Right under that, there should be some text saying "X items of data are available, covering Y participants, encoded using Data-Coding Z.") Click on the number "Z" to find out how the data is encoded. Put the value indicating "case" under `coding_binary_case` and the value indicating "control" under `coding_binary_control` in the Google Sheet.
  - According to the convention that we have, anything that is not case or control will automatically be encoded as missing. Put all such values in `coding_exclude`, delimited by commas.
- In the case of it being a quantitative trait, we need to look at how it is encoded. Some quantitative traits are encoded just as is, like [age started smoking in current smokers](http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=3436), which is encoded in the Biobank as years. If this is the case, leave all fields except name and field ID blank, but also include the special encodings into the Google sheet (for the link above, -1 means "Do not know" and -3 means "Prefer not to answer". These should be put into the `coding_exclude` category, semicolon-delimited). 
  - Some quantitative traits were measured using multiple choice questions. Imagine that you have a question like the following:
    
    Which of the following quantities typifies the quantity of sleep you get?
    - 0-2 hours
    - 2-4 hours
    - 4-6 hours
    - 6-8 hours
    - 8-10 hours
    - 10+ hours
    
    For this question, you should split your existing row into 6 rows. The naming convention for these kinds of rows is "Field Description (Parenthetical additional descriptor)" (e.g. "Average amount of sleep (0-2 hours)"). Additionally, there is an "ORDER" column in the sheet. For this kind of quantitative trait, the person doing this trait should ordinally assign values to this column based on quantity (1 for "0-2 hours", 2 for "2-4 hours", etc.).

Finally, you'll want an automated way of generating GBE IDs and field names for the phenotypes you want to annotate. After you have filled out the annotation table as above, the below macro can be helpful for doing just that. Paste this into the GBE ID column (column D) for the annotation at row X, replacing X with the row number, in order to automatically generate a GBE ID at row X:

```{excel}
=IF(SUM(GX,HX)>1, IF(ISBLANK(IX), CONCATENATE("BIN_FC",(JX*10000000)+FX), CONCATENATE("QT_FC",(IX*10000000)+FX)), IF(GX=1,CONCATENATE("INI",FX),CONCATENATE("BIN",FX)))
```

Likewise, you can automatically generate GBE-friendly field names at row X in column C:

```{excel}
=SUBSTITUTE(EX, " ", "_")
```

As in all Excel spreadsheets, you can strategically "drag down" this formula to make it populate across multiple rows if you're assigning IDs to multiple phenotypes during the session.
  
#### After the session - Compiling phenotype files
  
You are ready to now compile `.phe` files. This can be done in two ways. Export the updated Google Sheet as a `.tsv`, rename it as `ukb_YYYYMMDD.tsv`, and place it within [`tables/annotations`](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping/tables/annotations) (making sure to avoid namespace conflicts if applicable - it doesn't *have* to be the correct date, for example). Then, symlink [`tables/ukb_annotations.tsv`](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/tables/ukb_annotations.tsv) to this file (i.e., run `ln -sf tables/annotations/ukb_YYYYMMDD.tsv tables/ukb_annotations.tsv`). Then, continue with steps 8-10 as described [above](https://github.com/rivas-lab/ukbb-tools/tree/master/02_phenotyping#generating-and-updating-phenotypes-and-summary-statistics).
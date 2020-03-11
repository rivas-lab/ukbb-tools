# Phenotype files for snpnet-cox

This directory contains scripts and files for running snpnetcox. 

The location of the master snpnet-cox phenotype file for White British individuals:
`/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.cox.20200205.phe.zst`

The location of the master snpnet-cox phenotype file for all individuals:
`/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.cox.20200221.allinds.phe.zst`

## Compress files using zstd
```
$ ml load zstd
$ zstd -z file 
```

## Reshape OPCS4 codes
[Information on UKBB OPCS4 operative procedures](http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=41272).

The OPCS4 operative procedures are data-field 41272. Array indices run from 0-116.
There are matching data-fields with date of first operative procedure - data-field 41282 (0-116).
There is a total of 8125 different OPCD4 codes in a 4-character structure.
Script to reshape data format with columns FID and all the 8125 OPCS4 codes.

`opcs.formatting.R`

## OPCS4: Make individual phenotype files for snpnet-cox

Script

`opcs_snpnet_cox.R`

y - age at event or censoring/death

status - disease status: 0 is no event and 1 is event

inc - status on prevalent case = 0 and incident case = 1 (at UKBB baseline assessment)

## First Occurrence: Make individual phenotype files for snpnet-cox
[Documentation of UKBB first occurrence data-fields](http://biobank.ndph.ox.ac.uk/showcase/showcase/docs/first_occurrences_outcomes.pdf).
[UKBB first occurrence showcase data](http://biobank.ctsu.ox.ac.uk/crystal/search.cgi?wot=0&srch=first+occurrence&sta0=on&sta1=on&sta2=on&sta3=on&str0=on&str3=on&fit0=on&fit10=on&fit20=on&fit30=on&fvt11=on&fvt21=on&fvt22=on&fvt31=on&fvt41=on&fvt51=on&fvt61=on&fvt101=on).

The First Occurrence data-fields are composed of data from primary care, in-patient hospital, x and self-reported data. This information is all mapped to ICD-10 3-character structure with an event date. 


Script to make file with columns: FID, coxnet_y_(phenotype), coxnet_status_(phenotype) and coxnet_inc_(phenotype)
`first_occ_snpnet_cox.R`

## File on Sherlock
All individual phenotype files for First Occurrence (ICD-10 codes) and OPCS4 are stored in this directory on Sherlock:
`/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/opcs/phenotypefiles`


## Make new phenotype file based on multiple operational procedure codes (OPCS4)
Based on the files in the former section, we have a script in order to create new phenotype files by combining multiple OPCS4 codes. For some phenotypes we wish to combine several OPCS4 into one phenotype and use first event. 

To use disease diagnosis time to surgery as a proxy of disease progression.

Args[1] will subset to disease ICD-10 (status = 1)
Args[-1] are all the phenotypes you wish to combine as event.

Make a new phenotype file based on a disease (First Occurrence based on ICD10 codes) status = 1.
Script
`make_new_snpnet_cox_phe.R`

## Examples of running the script

```
$ Rscript make_new_snpnet_cox_phe.R f.130814.0.0 f.131306.0.0
```

The file generated includes all individuals who was diagnosed (status = 1) with f.130814.0.0 (disorders of lipoprotein metabolism and other lipidaemias|E78) and the event is f.131306.0.0 (chronic ischaemic heart disease|I25) 

Output file: `snpnet_cox_f.130814.0.0_f.131306.0.0.phe`

```
$ Rscript make_new_snpnet_cox_phe.R f.130708.0.0 X093 X094 X095
```

This file includes individuals having type 2 diabetes (non-insulin-dependent diabetes mellitus|E11) and the events are combined operational procedures for amputation of leg.

Output file: `snpnet_cox_f.130708.0.0_X093_X094_X095.phe`


## File that maps First Occurrence codes to phenotype name and ICD-10 code
Head of `mapfinal.txt`
``` 
filename|controls|cases|phenotypename|icd10
coxnet_status_f.131286.0.0.csv|253425|83727|essential (primary) hypertension|I10
coxnet_status_f.130814.0.0.csv|268675|68477|disorders of lipoprotein metabolism and other lipidaemias|E78
coxnet_status_f.131888.0.0.csv|279289|57863|other joint disorders, not elsewhere classified|M25
coxnet_status_f.131960.0.0.csv|282194|54958|other soft tissue disorders, not elsewhere classified|M79
```
Control and case numbers are for the White British population.


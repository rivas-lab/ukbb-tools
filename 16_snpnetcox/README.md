# Analyses snpnet-cox

This directory contains scripts and files for running snpnetcox. 

The location of the master cox phe file:
`/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.cox.20200205.phe`

## Reshape OPCS4 codes
[Information on UKBB OPCS4 operative procedures] (http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=41272)
The OPCS4 operative procedures are data-field 41272. Array indices run from 0-116.
There are matching data-fields with date of first operative procedure - data-field 41282 (0-116).
There is a total of 8125 different OPCD4 codes in a 4-character structure.
Script to reshape data format with columns FID and all the 8125 OPCS4 codes.
`opcs.formatting.R`

## Make individual phenotype files for snpnet-cox
[Documentation of UKBB first occurrence data-fields] (http://biobank.ndph.ox.ac.uk/showcase/showcase/docs/first_occurrences_outcomes.pdf)
[UKBB first occurrence showcase data] (http://biobank.ctsu.ox.ac.uk/crystal/search.cgi?wot=0&srch=first+occurrence&sta0=on&sta1=on&sta2=on&sta3=on&str0=on&str3=on&fit0=on&fit10=on&fit20=on&fit30=on&fvt11=on&fvt21=on&fvt22=on&fvt31=on&fvt41=on&fvt51=on&fvt61=on&fvt101=on)

Script to make file with columns: FID, coxnet_y_(phenotype), coxnet_status_(phenotype) and coxnet_inc_(phenotype)
`first_occ_snpnet_cox.R`

#!/bin/python

import pandas as pd

ftbd = pd.read_csv('../tables/field_table_basket_date.tsv', sep='\t')
dds = pd.read_csv('../tables/Data_Dictionary_Showcase.csv')
ann = pd.read_csv('../tables/ukb_annotations.tsv', sep='\t')

ftbd_set = set(ftbd['FieldID'])
dds_set = set(dds['FieldID'])
ann_set = set(ann['FieldID'])

all_fields = set(ftbd['FieldID']).union(set(dds['FieldID'])).union(set(ann['FieldID']))
print("There are " + str(len(dds_set)) + " fields in the Data Dictionary Showcase, " + str(len(ftbd_set)) + " fields in the ftbd tsv, and " + str(len(ann_set)) + " fields in our annotation file.")
print("There are " + str(len(all_fields)) + " fields total.")
print("We care about the " + str(len(ftbd_set)) + " fields in the ftbd tsv (all data we have gotten thus far).")

unannotated = ftbd_set - ann_set

print("Of these " +  str(len(ftbd_set)) + " fields in the ftbd tsv, " + str(len(unannotated)) + " are not in the annotation file.")
for field in unannotated:
    print(field)

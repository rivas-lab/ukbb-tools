#!/bin/python

import pandas as pd

ftbd = pd.read_csv('../tables/field_table_basket_date.tsv', sep='\t')
dds = pd.read_csv('../tables/Data_Dictionary_Showcase.csv')

diff = set(ftbd['FieldID']) - set(dds['FieldID'])

for field in diff:
    print(field)

print("")
print(str(len(diff)) + " total fields that we have that are not in Data Dictionary Showcase. Contact access@ukbiobank.ac.uk to get these into a new DDS before proceeding.")

#!/bin/python

import pandas as pd

rg = pd.read_table('/oak/stanford/groups/mrivas/ukbb24983/array-combined/ldsc/rg.metal.20201102-173634.tsv')
phen_info = pd.read_table('sumstat_paths.tsv')
phenos = list(phen_info['GBE_ID'])
rg = rg[(rg['p1'].isin(phenos)) & (rg['p2'].isin(phenos))]
rg = rg.fillna(0)
rg[['p1','p2','rg']].to_csv('biomarkers_array_rg.tsv', sep='\t', index=False)

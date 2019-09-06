#!/bin/python

import pandas as pd

tbl = pd.read_table('bbj_combined_sites.txt')

hello = tbl.groupby(['CHROM', 'POS', 'REF', 'ALT'], as_index=False).median()

hello.to_csv('bbj_combined_sites.txt', sep='\t', index=False)

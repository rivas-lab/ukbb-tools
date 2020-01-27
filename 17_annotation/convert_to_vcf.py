import sys

import pandas as pd
pd.options.mode.chained_assignment = None

fn = sys.argv[1]
out = sys.argv[2]
df = pd.read_table(fn)[['CHROM', 'POS', 'REF', 'ALT']]
df.columns = ['#CHROM', 'POS', 'REF', 'ALT']
df['FILTER'] = '.'
df['QUAL'] = '.'
df['INFO'] = '.'
df['ID'] = '.'
df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']]
df.to_csv(out, sep='\t', index=False)

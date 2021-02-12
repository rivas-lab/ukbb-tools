import pandas as pd

multipop = pd.read_table('mrp_rv_ma_exome_multipop_gbe.tsv')
singlepop = pd.read_table('mrp_rv_ma_exome_singlepop_gbe.tsv')

print("multipop")
print(len(multipop))
print("singlepop")
print(len(singlepop))
merged = pd.concat([multipop, singlepop], sort=False)
print(len(merged))
merged = merged.drop_duplicates()
print(len(merged))
merged.to_csv('mrp_rv_ma_exome_gbe.tsv', sep='\t', index=False)

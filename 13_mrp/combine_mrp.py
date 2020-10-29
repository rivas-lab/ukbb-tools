import pandas as pd

multipop = pd.read_table('mrp_rv_ma_array_multipop_gbe.txt')
singlepop = pd.read_table('mrp_rv_ma_array_singlepop_gbe.txt')

print(len(multipop))
print(len(singlepop))
merged = multipop.merge(singlepop, how='outer')
print(len(merged))
merged.to_csv('mrp_rv_ma_array_gbe.txt', sep='\t', index=False)

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

to_truncate = [col for col in merged.columns if (("log_10_BF" in col) or ("posterior_prob" in col))]
for col in to_truncate:
    merged[col] = merged[col].map("{0:.3g}".format)
merged.to_csv('mrp_rv_ma_exome_gbe.tsv', sep='\t', index=False)

import pandas as pd

df_37855 = pd.read_table('ukb_20191030.tsv', dtype=object)
df_37855 = df_37855.assign(Notes="")
df_37855 = df_37855.assign(Priority=2)
df_35059 = pd.read_table('ukb_20191101.2004890.35059.tsv', dtype=object)
df_35059 = df_35059.assign(Notes="")
df_35059 = df_35059.assign(Priority=2)
df_37338 = pd.read_table('ukb_20191101.2005223.37338.tsv', dtype=object)
df_37338 = df_37338.assign(Notes="")
df_37338 = df_37338.assign(Priority=2)

tables = ['ukb_20170727.tsv', 'ukb_20170818.tsv', 'ukb_20170827.tsv', 'ukb_20171015.tsv', 'ukb_20171110.tsv', 'ukb_20171113.tsv', 'ukb_20171211.tsv', 'ukb_20181109.tsv', 'ukb_20190327.tsv', 'ukb_20190406.tsv', 'ukb_20190409.tsv', 'ukb_20190412.tsv']

#tables = ['ukb_20191030.tsv']
dfs = [pd.read_table(table, dtype=object) for table in tables]
dfs = [df.assign(Priority=1) for df in dfs]
dfs.extend([df_37855, df_35059, df_37338])

big_df = pd.concat(dfs)
idx =  big_df.groupby(["FieldID"])["Priority"].transform(min) == big_df["Priority"]

df = big_df[idx]

df[["Annotator", "Annotation date", "Name", "GBE ID", "Field", "BasketID", "TableID", "FieldID", "QT_total_num", "BIN_total_num", "QT_index", "BIN_index", "coding_exclude", "coding_QT", "coding_binary_case", "coding_binary_control", "Participants", "Stability", "ValueType", "Units", "Strata", "Sexed", "Instances", "Array", "Coding", "Link", "Notes"]].to_csv('ukb_20191101_filled.tsv', sep='\t', index=False)

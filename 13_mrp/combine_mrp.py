import pandas as pd

multipop = pd.read_table('mrp_rv_ma_exome_multipop_gbe.tsv')
singlepop = pd.read_table('mrp_rv_ma_exome_singlepop_gbe.tsv')
biomarkers = list(pd.read_table('biomarkers_clusters.tsv')['GBE_ID'])

merged = pd.concat([multipop, singlepop], sort=False)
merged = merged.drop_duplicates()

def truncate(df, filename):
    to_truncate = [col for col in df.columns if (("log_10_BF" in col) or ("posterior_prob" in col))]
    for col in to_truncate:
        df[col] = df[col].map("{0:.3g}".format)
    df.to_csv(filename, sep='\t', index=False)

truncate(merged, 'mrp_rv_ma_exome_gbe.tsv')
biomarker_tbl = merged[merged['#GBE_ID'].isin(biomarkers)]
biomarker_tbl.to_csv('biomarkers_mrp_rv_ma_exome_gbe.tsv', sep='\t', index=False)

exome = pd.read_table('biomarkers_mrp_rv_exome_gbe.tsv')
truncate(exome, 'biomarkers_mrp_rv_exome_gbe.tsv')

array = pd.read_table('biomarkers_mrp_rv_ma_array_gbe.tsv')
truncate(array, 'biomarkers_mrp_rv_ma_array_gbe.tsv')

exome_var = pd.read_table('biomarkers_mrp_rv_exome_var_gbe.tsv')
truncate(exome_var, 'biomarkers_mrp_rv_exome_var_gbe.tsv')

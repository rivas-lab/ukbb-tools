import pandas as pd

e_multipop = pd.read_table('mrp_rv_ma_exome_multipop_gbe.tsv')
e_singlepop = pd.read_table('mrp_rv_ma_exome_singlepop_gbe.tsv')
e_rare_multipop = pd.read_table('mrp_rv_ma_exome_ultrarare_multipop_gbe.tsv')
e_rare_singlepop = pd.read_table('mrp_rv_ma_exome_ultrarare_singlepop_gbe.tsv')
a_multipop = pd.read_table('mrp_rv_ma_array_multipop_gbe.tsv')
a_singlepop = pd.read_table('mrp_rv_ma_array_singlepop_gbe.tsv')
biomarkers = list(pd.read_table('biomarkers_clusters.tsv')['GBE_ID'])

e_merged = pd.concat([e_multipop, e_singlepop], sort=False)
e_merged = e_merged.drop_duplicates()
e_rare_merged = pd.concat([e_rare_multipop, e_rare_singlepop], sort=False)
e_rare_merged = e_rare_merged.drop_duplicates()
a_merged = pd.concat([a_multipop, a_singlepop], sort=False)
a_merged = a_merged.drop_duplicates()

def truncate(df, filename):
    to_truncate = [col for col in df.columns if (("log_10_BF" in col) or ("posterior_prob" in col))]
    for col in to_truncate:
        df[col] = df[col].map("{0:.3g}".format)
    df.to_csv(filename, sep='\t', index=False)

truncate(e_merged, 'mrp_rv_ma_exome_gbe.tsv')
truncate(e_rare_merged, 'mrp_rv_ma_exome_ultrarare_gbe.tsv')
truncate(a_merged, 'mrp_rv_ma_array_gbe.tsv')

e_biomarker_tbl = e_merged[e_merged['#GBE_ID'].isin(biomarkers)]
e_biomarker_tbl.to_csv('biomarkers_mrp_rv_ma_exome_gbe.tsv', sep='\t', index=False)
e_biomarker_rare_tbl = e_rare_merged[e_rare_merged['#GBE_ID'].isin(biomarkers)]
e_biomarker_rare_tbl.to_csv('biomarkers_mrp_rv_ma_exome_ultrarare_gbe.tsv', sep='\t', index=False)
a_biomarker_tbl = a_merged[a_merged['#GBE_ID'].isin(biomarkers)]
a_biomarker_tbl.to_csv('biomarkers_mrp_rv_ma_array_gbe.tsv', sep='\t', index=False)

exome = pd.read_table('biomarkers_mrp_rv_exome_gbe.tsv')
truncate(exome, 'biomarkers_mrp_rv_exome_gbe.tsv')

exome_var = pd.read_table('biomarkers_mrp_rv_exome_var_gbe.tsv')
truncate(exome_var, 'biomarkers_mrp_rv_exome_var_gbe.tsv')

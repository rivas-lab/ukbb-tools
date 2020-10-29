import pandas as pd
import sys
import numpy as np

gbe = pd.read_table('mrp_rv_ma_array_gbe.txt')
phen_info = pd.read_table('../05_gbe/array-combined/phenotype_info.tsv')[['#GBE_ID','GBE_NAME']]
gbe = gbe.merge(phen_info, on='#GBE_ID')

genes = np.unique(gbe['gene'])

for gene in genes:
    print(gene)
    subset_df = gbe[gbe['gene'] == gene][['#GBE_ID', 'GBE_NAME', 'log_10_BF_study_similar_var_similar_sigma_m_var_ptv', 'log_10_BF_study_similar_var_independent_sigma_m_var_pav', 'log_10_BF_study_similar_var_similar_sigma_m_var_pav']]
    if sys.argv[1] == "bf_filter":
        subset_df = subset_df[(subset_df['log_10_BF_study_similar_var_similar_sigma_m_var_ptv'] >= .5) | (subset_df['log_10_BF_study_similar_var_independent_sigma_m_var_pav'] >= .5) | (subset_df['log_10_BF_study_similar_var_similar_sigma_m_var_pav'] >= .5)]
        subset_df.to_csv('reformat_mrp/' + gene + '_bf_filter.tsv', index=False)
    else:
        subset_df.to_csv('reformat_mrp/' + gene + '.tsv', index=False)

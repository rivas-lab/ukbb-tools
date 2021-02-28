#!/bin/python

import pandas as pd
import numpy as np
import sys

results = pd.read_table(sys.argv[1])
biomarkers = pd.read_table('biomarkers.tsv')
results = results[results['GBE_short_name'].isin(list(biomarkers['#Phenotype']))]

to_replace = ['log_10_BF_study_similar_var_independent_sigma_m_mpc_pli_pav', 'posterior_prob_w_prior_odds_0.0005_study_similar_var_independent_sigma_m_mpc_pli_pav', 'log_10_BF_study_similar_var_independent_sigma_m_mpc_pli_ptv', 'posterior_prob_w_prior_odds_0.0005_study_similar_var_independent_sigma_m_mpc_pli_ptv']
replacements = ['log_10_BF_study_independent_var_independent_sigma_m_mpc_pli_pav', 'posterior_prob_w_prior_odds_0.0005_study_independent_var_independent_sigma_m_mpc_pli_pav', 'log_10_BF_study_independent_var_independent_sigma_m_mpc_pli_ptv', 'posterior_prob_w_prior_odds_0.0005_study_independent_var_independent_sigma_m_mpc_pli_ptv']
for col1, col2 in zip(to_replace, replacements):
    results[col1] = np.where(results['num_pops'] == 1, results[col2], results[col1])

pav = results[results['log_10_BF_study_similar_var_independent_sigma_m_mpc_pli_pav'] >= 5]
ptv = results[results['log_10_BF_study_similar_var_independent_sigma_m_mpc_pli_ptv'] >= 5]

"""
print("PAV num associations")
print(len(pav))
print("num traits")
print(len(pav['#GBE_ID'].unique()))
print("num genes")
print(len(pav['gene'].unique()))
print("PTV num associations")
print(len(ptv))
print("num traits")
print(len(ptv['#GBE_ID'].unique()))
print("num genes")
print(len(ptv['gene'].unique()))


"""
gene = pd.read_table('ref/gene_chr_pos.tsv')
pav = pav.merge(gene, how='left')
ptv = ptv.merge(gene, how='left')

to_include = pd.read_table('pheno_fields.tsv')

colors = ['#fc0303', '#fc7b03', '#fce703', '#67fc03', '#03fcfc', '#037bfc', '#030bfc', '#8403fc', '#eb03fc', '#8c0000', '#8c4800', '#8c8300', '#368c00', '#008a8c', '#002c8c', '#23008c', '#67008c', '#8c006e', '#ff8787', '#ffc387', '#fff587', '#c1ff87', '#87ffa3', '#87fffb', '#87c7ff', '#8789ff', '#c787ff', '#f787ff', '#ff87cb']

for name, df in zip(['pav', 'ptv'], [pav, ptv]):
    df = df[df['#GBE_ID'].isin(list(to_include['title']))]
    df = df.assign(MARKER=(df['CHROM'].astype(str) + '_' + df['POS'].astype(str) + '_' + df['gene']))
    df = df.sort_values(['CHROM', 'POS'])
    df = df.assign(LOCUS_ID=df['MARKER'].astype('category').cat.codes)
    df['LOCUS_ID'] = df['LOCUS_ID'].astype(int) + 1
    df = df[['LOCUS_ID', 'GBE_short_name', 'CHROM', 'POS', 'MARKER', 'gene']]
    df.columns = ['LOCUS_ID', 'TRAIT', 'CHR', 'BP', 'MARKER', 'GENE']
    color_df = biomarkers[['CATEGORY', '#Phenotype', 'COLOR']].drop_duplicates()
    color_df.columns = ['CATEGORY', 'TRAIT', 'COLOR']
    df = df.merge(color_df[['CATEGORY','TRAIT']])
    df = df[['LOCUS_ID', 'CATEGORY', 'TRAIT', 'CHR', 'BP', 'MARKER', 'GENE']]
    df.to_csv(name + '_fuji.tsv', sep='\t', index=False)
    #categories = df['CATEGORY'].unique()
    #num_categories = len(categories)
    #colors = colors[:num_categories]
    #color_dict = dict(zip(categories, colors))
    #color_df['COLOR'] = color_df['CATEGORY'].map(color_dict)
    #color_df = color_df.sort_values(['CATEGORY', 'TRAIT'])
    color_df.to_csv(name + '_color.tsv', sep='\t', index=False)

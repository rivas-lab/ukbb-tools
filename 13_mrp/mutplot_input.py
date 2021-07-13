import pandas as pd

exm = pd.read_table('ref/ukb_exm_oqfe-consequence_wb_maf_gene_ld_indep_mpc_pli.tsv')
#exm = exm[exm['MPC'].notna()]
alpl = exm[exm['gene_symbol'] == 'ALPL']

for df, filename in zip([alpl], ['alpl.tsv']):
    df['protein_change'] = df.HGVSp.str.split('p.').str[1]
    df['aa'] = df.protein_change.str.extract('(\d+)')
    df = df[df['aa'].notna()]
    df['aa'] = df['aa'].astype(int)
    df = df[['V', 'gene_symbol', 'protein_change', 'aa', 'most_severe_consequence', 'MPC']]
    df = df[df['MPC'].notna()]
    df['most_severe_consequence'] = df['most_severe_consequence'].str.replace('_', ' ')
    df.to_csv(filename, sep='\t', index=False)

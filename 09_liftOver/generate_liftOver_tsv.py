import pandas as pd

bim_df = pd.read_csv(
    'ukb24983_cal_cALL_v2.bim',
    sep='\t', names=['#CHROM', 'ID', 'dist', 'hg37_pos', 'hg37_A1', 'hg37_A2']
)

liftOver_df = pd.read_csv(
    'hg38', sep='\t', usecols=[1,3], names=['hg38_pos', 'ID']
)

df = bim_df.merge(
    liftOver_df, how='left', on='ID'
).fillna(value={'hg38_pos':0})

df['hg38_pos'] = df['hg38_pos'].map(lambda x: int(x))

df[['#CHROM', 'ID', 'hg37_pos', 'hg38_pos', 'hg37_A1', 'hg37_A2']].to_csv(
    'ukb24983_cal_cALL_v2.liftOver.tsv',
    sep='\t', index=False
)


import pandas as pd

met = pd.read_table('mrp_rv_metabolomics.tsv')

def truncate(df, filename):
    to_truncate = [col for col in df.columns if (("log_10_BF" in col) or ("posterior_prob" in col))]
    for col in to_truncate:
        df[col] = df[col].map("{0:.3g}".format)
    df.to_csv(filename, sep='\t', index=False)

truncate(met, 'mrp_rv_metabolomics.tsv')

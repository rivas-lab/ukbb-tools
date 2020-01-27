import sys
import pandas as pd
import gzip

tsv = sys.argv[1]
vcf = sys.argv[2]
maf = sys.argv[3]
out = sys.argv[4]

tsv_df = pd.read_table(tsv)

i = 0

with gzip.open(vcf,'rt') as f:
    for line in f:
        if line.strip().split()[0] == "#CHROM":
            break
        else:
            i += 1

vcf_df = pd.read_table(vcf, skiprows=i)[['#CHROM', 'POS', 'REF', 'ALT', 'INFO']]
vcf_df.columns = ['CHROM', 'POS', 'REF', 'ALT', 'consequence_field']

maf_df = pd.read_table(maf)[['CHROM', 'POS', 'REF', 'ALT', 'MAF']]

out_df = tsv_df.merge(vcf_df)
out_df = out_df.merge(maf_df)

out_df = out_df.replace(pd.np.nan, "NA", regex=True)
out_df.to_csv(out, sep='\t', index=False)

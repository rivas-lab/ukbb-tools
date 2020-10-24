import pandas as pd
import re
import glob
import os

info = pd.read_table('exome_phenotype_info.tsv')
shortnames = pd.read_table('icdinfo.shortnames.exome.tsv')
info = info.merge(shortnames)

"""
for pop, col in zip(['white_british', 'non_british_white', 's_asian', 'e_asian', 'african', 'related', 'others'], ['N_GBE', 'N_NBW', 'N_SAS', 'N_EAS', 'N_AFR', 'N_SMR', 'N_OTH']):
    outfile = 'icdinfo.exome.' + pop +  '.tsv'
    data = []
    for i, row in info.iterrows():
        category = re.sub('\d', '', row['#GBE_ID'])
        data.append([category, row['#GBE_ID'], row[col], row['GBE_NAME'], row['GBE_short_name'], row['GBE_short_name_len']])
    to_write = pd.DataFrame(data, columns=['GBE_category', 'GBE_ID', 'N', 'GBE_NAME', 'GBE_short_name', 'GBE_short_name_len'])
    to_write.to_csv(outfile, sep='\t', index=False)
"""

exome_metals = glob.glob('/oak/stanford/groups/mrivas/ukbb24983/exome/metal/20201014/*txt')
gbe_ids = [os.path.basename(filepath).split(".")[0] for filepath in exome_metals]
num_dict = {'white_british': 'N_GBE',
            'non_british_white': 'N_NBW',
            's_asian': 'N_SAS',
            'e_asian': 'N_EAS',
            'african': 'N_AFR',
            'related': 'N_SMR',
            'others': 'N_OTH'}

data = []
for gbe_id, exome_metal in zip(gbe_ids, exome_metals):
    with open(exome_metal) as f:
        content = f.readlines()
    content = [x.strip() for x in content if ("oak" in x) and ("metal" not in x)]
    pops = [x.split("/")[-2] for x in content]
    cols = [num_dict[pop] for pop in pops]
    row = info[info['#GBE_ID'] == gbe_id]
    nums_metal = [row[col].values[0] for col in cols]
    num_metal = sum(nums_metal)
    data.append([re.sub('\d', '', gbe_id), gbe_id, num_metal, row['GBE_NAME'].values[0], row['GBE_short_name'].values[0], row['GBE_short_name_len'].values[0]])

to_write = pd.DataFrame(data, columns=['GBE_category', 'GBE_ID', 'N', 'GBE_NAME', 'GBE_short_name', 'GBE_short_name_len'])
to_write.to_csv('icdinfo.exome.metal.tsv', sep='\t', index=False)

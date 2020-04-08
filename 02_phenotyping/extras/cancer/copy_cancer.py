import pandas as pd
import numpy as np
import os
from shutil import copyfile

cancer = pd.read_table('cancer_gbe_map.tsv', header=None, names=['GBE_ID', 'phenotype'])
cancer['ID'] = cancer['GBE_ID'].str.replace('cancer','')

shitty_prefix = '/oak/stanford/groups/mrivas/users/mrivas/repos/ukbb-phenotyping/cancer'
good_prefix = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/cancer/v2_2020/phe/'

shitty_files = []

for i, row in cancer.iterrows():
    shitty_file = os.path.join(shitty_prefix, row['ID'] + '.phe')
    good_file = os.path.join(good_prefix, row['GBE_ID'] + '.phe')
    if not os.path.exists(shitty_file):
        print("BAD FILE: " + shitty_file)
        shitty_files.append(row['shitty_hc'])    
    else:
        copyfile(shitty_file, good_file)
    print(shitty_file, good_file)

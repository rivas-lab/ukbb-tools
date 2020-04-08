import pandas as pd
import numpy as np
import os
from shutil import copyfile

dst_prefix = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/highconfidenceqc/v2_2020/phe/'

tte = pd.read_table('TTE_HC.tsv')
ad = pd.read_table('AD_HC.tsv')

def copy_files(df):
    for i, row in df.iterrows():
        src_file = row['Source']
        dst_file = os.path.join(dst_prefix, row['GBE ID'] + '.phe')
        if not os.path.exists(src_file):
            print("BAD FILE: " + src_file)
        else:
            copyfile(src_file, dst_file)
        print(src_file, dst_file)

copy_files(ad)
copy_files(tte)

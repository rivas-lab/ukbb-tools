import os
import glob

qts = glob.glob('/oak/stanford/groups/mrivas/bbj/raw_sumstats/qt')
bins = glob.glob('/oak/stanford/groups/mrivas/bbj/raw_sumstats/bin')
eqtls = glob.glob('/oak/stanford/groups/mrivas/bbj/raw_sumstats/eqtl')

with open('bbj_name_map.tsv', 'r') as map_file:
    prefix = '/oak/stanford/groups/mrivas/bbj/raw_sumstats/'
    for line in map_file:
        if (not 'Folder' in line) and (not 'eqtl' in line):
            split_line = line.split()
            folder_name = split_line[0]
            curr_name = split_line[1] + '.txt.gz'
            target_name = os.path.join(prefix, folder_name, split_line[2])
            print(os.path.join(prefix, folder_name, curr_name), target_name)
            os.rename(os.path.join(prefix, folder_name, curr_name), target_name)

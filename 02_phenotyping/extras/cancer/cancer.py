import os
from annotate_phe import make_phe_info
import glob
import pandas as pd

if __name__ == "__main__":
    cancer_phenos = glob.glob('/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/cancer/phe/*.phe')
    with open('cancer_gbe_map.tsv', 'r') as map_file:
        gbe_id_to_name = {line.split()[0]:line.split()[1] for line in map_file}
    gbe_ids = [os.path.splitext(os.path.basename(path))[0] for path in cancer_phenos]
    cancer_names = [gbe_id_to_name[gbe_id] for gbe_id in gbe_ids]
    make_phe_info(in_phe   = cancer_phenos,
                  out_path = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/cancer/info/',
                  name     = cancer_names,
                  field    = 'NA',
                  table    = 'NA',
                  basket   = 'NA',
                  app_id   = '24983',
                  source   = 'NA',
                  one_file = False)

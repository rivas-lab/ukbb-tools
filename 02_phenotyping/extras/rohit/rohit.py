from annotate_phe import make_phe_info
import os
import glob
import pandas as pd

if __name__ == "__main__":
    rh_phenos = glob.glob('/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/rohit/phe/*.phe')
    with open('rohit_gbe_map.tsv', 'r') as map_file:
        gbe_id_to_name = {line.split()[0]:line.split()[1] for line in map_file}
    gbe_ids = [os.path.splitext(os.path.basename(path))[0] for path in rh_phenos]
    rh_names = [gbe_id_to_name[gbe_id] for gbe_id in gbe_ids]
    make_phe_info(in_phe   = rh_phenos,
                  out_path = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/cancer/info/',
                  name     = rh_names,
                  field    = 'NA',
                  table    = 'NA',
                  basket   = 'NA',
                  app_id   = '24983',
                  source   = 'rohit',
                  one_file = False)

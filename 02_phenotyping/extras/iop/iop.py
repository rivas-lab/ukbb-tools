import os
from annotate_phe import make_phe_info
import glob

if __name__ == "__main__":
    iop_phenos = glob.glob('/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/iop/phe/*.phe')
    with open('iop_gbe_map.tsv', 'r') as map_file:
        gbe_id_to_name = {line.split()[0]:line.split()[1] for line in map_file} 
    gbe_ids = [os.path.splitext(os.path.basename(path))[0] for path in iop_phenos]
    iop_names = [gbe_id_to_name[gbe_id] for gbe_id in gbe_ids]
    make_phe_info(in_phe   = iop_phenos,
                  out_path = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/iop/info/',
                  name     = iop_names,
                  field    = '5254/5255',
                  table    = 'NA',
                  basket   = 'NA',
                  app_id   = '24983',
                  source   = 'ytanigaw',
                  one_file = False)

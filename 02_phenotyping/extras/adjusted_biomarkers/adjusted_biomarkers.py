import os
from annotate_phe import make_phe_info
import glob
import pandas as pd

if __name__ == "__main__":
    ab_phenos = glob.glob('/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/adjusted_biomarkers/phe/*.phe')
    gbe_ids = [os.path.splitext(os.path.basename(filepath))[0] for filepath in ab_phenos]
    with open('adjusted_biomarker_gbe_map.tsv', 'r') as gbe_map:
        gbe_to_ab = {line.split()[0]:line.split()[1] for line in gbe_map} 
    ab_names = [gbe_to_ab[gbe_id] for gbe_id in gbe_ids]
    make_phe_info(in_phe   = ab_phenos,
                  out_path = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/adjusted_biomarkers/info/',
                  name     = ab_names,
                  field    = 'NA',
                  table    = 'NA',
                  basket   = 'NA',
                  app_id   = '24983',
                  source   = 'nasa',
                  one_file = False)

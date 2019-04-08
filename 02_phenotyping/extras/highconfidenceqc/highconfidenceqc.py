from annotate_phe import make_phe_info
import os
import glob

if __name__ == "__main__":
    hc_phenos  = glob.glob('/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/highconfidenceqc/phe/*.phe')
    with open('highconfidenceqc_gbe_map.tsv', 'r') as map_file:
        gbe_id_to_name = {line.split()[0]:line.split()[1] for line in map_file}
    gbe_ids = [os.path.splitext(os.path.basename(path))[0] for path in hc_phenos]
    hc_names = [gbe_id_to_name[gbe_id] for gbe_id in gbe_ids]
    make_phe_info(in_phe   = hc_paths,
                  out_path = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/highconfidenceqc/info/',
                  name     = hc_names,
                  field    = 'NA',
                  table    = 'NA',
                  basket   = 'NA',
                  app_id   = '24983',
                  source   = 'cdeboever',
                  one_file = False)

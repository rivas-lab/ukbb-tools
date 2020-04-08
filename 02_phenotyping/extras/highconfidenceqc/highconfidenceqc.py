from annotate_phe import make_phe_info
import os
import glob

if __name__ == "__main__":
    hc_phenos  = glob.glob('/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/highconfidenceqc/current/phe/*.phe')
    hc_phenos = [os.path.realpath(hc_pheno) for hc_pheno in hc_phenos]
    with open('highconfidenceqc_gbe_map.tsv', 'r') as map_file:
        gbe_id_to_name = {line.split()[0]:line.split()[1] for line in map_file}
    gbe_ids = [os.path.splitext(os.path.basename(path))[0] for path in hc_phenos]
    hc_names = [gbe_id_to_name[gbe_id] for gbe_id in gbe_ids]
    make_phe_info(in_phe   = hc_phenos,
                  out_path = os.path.realpath('/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/highconfidenceqc/current/info/'),
                  name     = hc_names,
                  field    = 'NA',
                  table    = 'NA',
                  basket   = 'NA',
                  app_id   = '24983',
                  source   = 'cdeboever/justesen',
                  one_file = False)

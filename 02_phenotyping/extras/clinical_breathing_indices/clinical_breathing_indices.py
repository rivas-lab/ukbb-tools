from annotate_phe import make_phe_info
import os
import glob

if __name__ == "__main__":
    cbi_phenos = glob.glob('/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/clinical_breathing_indices/phe/*.phe')
    with open('clinical_breathing_indices.tsv', 'r') as map_file:
        gbe_id_to_name = {line.split()[0]:line.split()[1] for line in map_file}
    gbe_ids = [os.path.splitext(os.path.basename(path))[0] for path in cbi_phenos]
    cbi_names = [gbe_id_to_name[gbe_id] for gbe_id in gbe_ids]
    make_phe_info(in_phe   = cbi_phenos,
                  out_path = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/clinical_breathing_indices/info/',
                  name     = cbi_names,
                  field    = 'NA',
                  table    = 'NA',
                  basket   = 'NA',
                  app_id   = '24983',
                  source   = 'magu',
                  one_file = False)

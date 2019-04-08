from annotate_phe import make_phe_info
import os
import glob

if __name__ == "__main__":
    with open('./clinical_breathing_indices_icdinfo.txt', 'r') as icd:
        icdinfo = {line.split()[0]:line.split()[2] for line in icd}
    cbi_phenos = glob.glob('/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/clinical_breathing_indices/phe/*.phe')
    cbi_names  = list(map(lambda c: icdinfo[c] if c in icdinfo else '',
                         map(lambda path: os.path.splitext(os.path.basename(path))[0], cbi_phenos)))
    make_phe_info(in_phe   = cbi_phenos,
                  out_path = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/clinical_breathing_indices/info/',
                  name     = cbi_names,
                  field    = 'NA',
                  table    = 'NA',
                  basket   = 'NA',
                  app_id   = '24983',
                  source   = 'magu',
                  one_file = False)

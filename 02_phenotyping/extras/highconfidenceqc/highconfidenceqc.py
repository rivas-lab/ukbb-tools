from annotate_phe import make_phe_info
import os
import glob

if __name__ == "__main__":
    with open('../../../../wiki/ukbb/icdinfo/icdinfo.txt', 'r') as icd:
        icdinfo = {line.split()[0]:line.split()[2] for line in icd}
    hc_paths  = glob.glob('/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/highconfidenceqc/phe/*.phe')
    hc_names  = list(map(lambda hc: icdinfo[hc] if hc in icdinfo else '',
                         map(lambda path: os.path.splitext(os.path.basename(path))[0], hc_paths)))
    make_phe_info(in_phe   = hc_paths,
                  out_path = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/highconfidenceqc/info/',
                  name     = hc_names,
                  field    = 'NA',
                  table    = 'NA',
                  basket   = 'NA',
                  app_id   = '24983',
                  source   = 'cdeboever',
                  one_file = False)

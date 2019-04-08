from annotate_phe import make_phe_info
import os
import glob

if __name__ == "__main__":
    fh_phenos = glob.glob('/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/family_history/phe/*.phe')
    fh_map    = {'FH'+line.split()[0][:4]:line.rstrip().split(None,1)[1].replace(' ','_') for line in open('/oak/stanford/groups/mrivas/private_data/ukbb/16698/phenotypedata/familyHistory2/mapphe.txt', 'r')}
    fh_names  = [fh_map[phe[:6]] if phe[:6] in fh_map else '' for phe in map(os.path.basename,fh_phenos)]
    make_phe_info(in_phe   = fh_phenos,
                  out_path = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/family_history/info/',
                  name     = fh_names,
                  field    = 'NA',
                  table    = 'NA',
                  basket   = 'NA',
                  app_id   = '24983',
                  source   = 'rohit',
                  one_file = False)

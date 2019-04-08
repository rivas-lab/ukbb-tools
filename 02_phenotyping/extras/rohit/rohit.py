from annotate_phe import make_phe_info
import os
import glob

if __name__ == "__main__":
    rh_phenos = glob.glob('/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/rohit/phe/*.phe')
    rh_numbers = ['RH' + line.rstrip().split()[0] for line in open('rohit_map.txt', 'r')]
    rh_names  = [('_').join(line.rstrip().split()[1:]) for line in open('rohit_map.txt', 'r')]
    rh_dict = dict(zip(rh_numbers, rh_names))    
    rh_ordered_phenos = [os.path.splitext(os.path.basename(path))[0] for path in rh_phenos]
    rh_ordered_names = [rh_dict[pheno] for pheno in rh_ordered_phenos]
    make_phe_info(in_phe   = rh_phenos,
                  out_path = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/rohit/info/',
                  name     = rh_ordered_names,
                  field    = 'NA',
                  table    = 'NA',
                  basket   = 'NA',
                  app_id   = '24983',
                  source   = 'rohit',
                  one_file = False)

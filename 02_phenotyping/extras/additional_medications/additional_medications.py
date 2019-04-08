from annotate_phe import make_phe_info
import os
import glob

if __name__ == "__main__":
    med_phenos = glob.glob('/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/additional_medications/phe/*.phe')
    med_names = [line.rstrip().split()[0] for line in open('med_names.txt', 'r')]
    med_numbers  = [line.rstrip().split()[1] for line in open('med_names.txt', 'r')]
    med_dict = dict(zip(med_numbers, med_names))
    med_ordered_phenos = [os.path.splitext(os.path.basename(path))[0] for path in med_phenos]
    med_ordered_names = [med_dict[pheno] for pheno in med_ordered_phenos]
    make_phe_info(in_phe   = med_phenos,
                  out_path = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/additional_medications/info/',
                  name     = med_ordered_names,
                  field    = 'NA',
                  table    = 'NA',
                  basket   = 'NA',
                  app_id   = '24983',
                  source   = 'alavertu',
                  one_file = False)

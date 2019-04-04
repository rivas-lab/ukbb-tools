from annotate_phe import make_phe_info

if __name__ == "__main__":
    import glob
    with open('../../../../wiki/ukbb/icdinfo/icdinfo.txt', 'r') as icd:
        icdinfo = {line.split()[0]:line.split()[2] for line in icd}
    cancer_phenos = glob.glob('/oak/stanford/groups/mrivas/private_data/ukbb/16698/phenotypedata/cancer3/*.phe')
    cancer_names  = list(map(lambda c: icdinfo['cancer'+c] if 'cancer'+c in icdinfo else '',
                         map(lambda path: os.path.splitext(os.path.basename(path))[0], cancer_phenos)))
    make_phe_info(in_phe   = cancer_phenos,
                  out_path = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/cancer/info/',
                  name     = cancer_names,
                  field    = 'NA',
                  table    = 'NA',
                  basket   = 'NA',
                  app_id   = '16698',
                  source   = 'NA',
                  one_file = False)

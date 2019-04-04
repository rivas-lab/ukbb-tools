from annotate_phe import make_phe_info

if __name__ == "__main__":
    import glob
    with open('../../../../wiki/ukbb/icdinfo/icdinfo.txt', 'r') as icd:
        icdinfo = {line.split()[0]:line.split()[2] for line in icd}
    fh_phenos = glob.glob('/oak/stanford/groups/mrivas/private_data/ukbb/16698/phenotypedata/familyHistory2/FH*.phe')
    fh_map    = {'FH'+line.split()[0][:4]:line.rstrip().split(None,1)[1].replace(' ','_') for line in open('/oak/stanford/groups/mrivas/private_data/ukbb/16698/phenotypedata/familyHistory2/mapphe.txt', 'r')}
    fh_names  = [fh_map[phe[:6]] if phe[:6] in fh_map else '' for phe in map(os.path.basename,fh_phenos)]
    make_phe_info(in_phe   = fh_phenos,
                  out_path = '/oak/stanford/groups/mrivas/dev-ukbb-tools/phenotypes/extras/family_history/info/',
                  name     = fh_names,
                  field    = 'NA',
                  table    = 'NA',
                  basket   = 'NA',
                  app_id   = '16698',
                  source   = 'rohit',
                  one_file = False)

from annotate_phe import make_phe_info

if __name__ == "__main__":
    import glob
    with open('../../../../wiki/ukbb/icdinfo/icdinfo.txt', 'r') as icd:
        icdinfo = {line.split()[0]:line.split()[2] for line in icd}
    hc_paths  = glob.glob('/oak/stanford/groups/mrivas/private_data/ukbb/16698/phenotypedata/highconfidenceqc/HC*.phe')
    hc_names  = list(map(lambda hc: icdinfo[hc] if hc in icdinfo else '',
    make_phe_info(in_phe   = hc_paths,
                  out_path = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/highconfidenceqc/info/',
                  name     = hc_names,
                  field    = 'NA',
                  table    = 'NA',
                  basket   = 'NA',
                  app_id   = '16698',
                  source   = 'cdeboever',
                  one_file = False)

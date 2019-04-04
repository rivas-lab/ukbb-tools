from annotate_phe import make_phe_info

if __name__ == "__main__":
    import glob
    with open('../../../../wiki/ukbb/icdinfo/icdinfo.txt', 'r') as icd:
        icdinfo = {line.split()[0]:line.split()[2] for line in icd}
    rh_phenos = ['/oak/stanford/groups/mrivas/private_data/ukbb/16698/phenotypedata/rohit/RH{}.phe'.format(i) for i in range(161)]
    rh_names  = [line.rstrip().split()[1] for line in open('/oak/stanford/groups/mrivas/private_data/ukbb/16698/phenotypedata/rohit/rohit_map.txt', 'r')]
    make_phe_info(in_phe   = rh_phenos,
                  out_path = '/oak/stanford/groups/mrivas/dev-ukbb-tools/phenotypes/extras/rohit/info/',
                  name     = rh_names,
                  field    = 'NA',
                  table    = 'NA',
                  basket   = 'NA',
                  app_id   = '16698',
                  source   = 'rohit',
                  one_file = False)

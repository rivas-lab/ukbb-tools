#!/bin/python
import os
from datetime import date

def make_phe_info(in_phe, out_path, name, field, table, basket, app_id, source, one_file=False):
    # reference 
    pop_dir = '/oak/stanford/groups/mrivas/private_data/ukbb/' + app_id + '/sqc/population_stratification/'
    get_iids = lambda f: set([line.split()[0] for line in open(f, 'r')])
    pops = {pop:get_iids(pop_dir + 'ukb'+ app_id + '_' + pop + '.phe')  
                for pop in ('white_british','african','e_asian','s_asian')}
    # this is what gets written to file, one item per line below 
    header = "#GBE_ID GBE_NAME FIELD TABLE BASKET APP_ID N N_GBE N_AFR N_EAS N_SAS SOURCE DATE PATH".replace(" ", "\t")
    # in_phe and name are lists of equal sizes -- everything else is constant strings
    if one_file:
        o = open(out_path, 'w')
        o.write(header + '\n')
    # iterate over paths to files and their names (not* GBE IDs)
    for phe_path, gbe_name in zip(in_phe, name):
        # this is a lazy input filtering mechanism from main()
        if not gbe_name: continue
        # load phenotype data 
        phe = dict([line.split()[1:] for line in open(phe_path, 'r') if 'IID' not in line])
        gbe_id = os.path.splitext(os.path.basename(phe_path))[0]
        # more filtering
        if 'cancer3' in phe_path:
            gbe_id = 'cancer' + gbe_id
        # determine which values count towards N (assuming plink formatted file)
        # where case/control is 2/1 and missing is -9
        if all([i in ['1','2','-9'] for i in phe.values()]): # binary
            count = lambda x: x == '2' 
        else: # quantitative
            count = lambda x: float(x) != -9
        # get info — each line (approximately) corresponds to one item in the header above
        phe_info = '\t'.join([gbe_id,
                              gbe_name.replace(' ','_'),
                              field,
                              table,
                              basket,
                              app_id,
                              str(sum(map(count,phe.values())))] +
                             [str(sum(map(lambda x:count(x[1]),
                                          filter(lambda x:x[0] in pops[pop], phe.items()))))
                              for pop in ('white_british','african','e_asian','s_asian')] +                   
                             [source,
                              str(date.today()),
                              os.path.abspath(phe_path)]) + '\n'
        # manage file writing options
        if not one_file:
            o = open(os.path.join(out_path, gbe_id + '.info'), 'w')
            o.write(header + '\n')
        o.write(phe_info)
        if not one_file:
            o.close()
        # progress bar, effectively — this should be kept in case it's used somewhere else
        print(phe_info[:-1])
    if one_file:
        o.close()
    return


if __name__ == "__main__":
    import glob
    # -1: load old icdinfo for reference
    with open('../../../wiki/ukbb/icdinfo/icdinfo.txt', 'r') as icd:
        icdinfo = {line.split()[0]:line.split()[2] for line in icd}
    # 0. write out HC info files to master table
    hc_paths  = glob.glob('/oak/stanford/groups/mrivas/private_data/ukbb/16698/phenotypedata/highconfidenceqc/HC*.phe')
    hc_names  = list(map(lambda hc: icdinfo[hc] if hc in icdinfo else '', 
                         map(lambda path: os.path.splitext(os.path.basename(path))[0], hc_paths)))
    make_phe_info(in_phe   = hc_paths, 
                  out_path = '/oak/stanford/groups/mrivas/dev-ukbb-tools/phenotypes/extra_info/hc.tsv.info',
                  name     = hc_names,
                  field    = 'NA', 
                  table    = 'NA', 
                  basket   = 'NA',
                  app_id   = '16698', 
                  source   = 'cdeboever', 
                  one_file = True)
    # 1. write out cancer info files to master table
    cancer_phenos = glob.glob('/oak/stanford/groups/mrivas/private_data/ukbb/16698/phenotypedata/cancer3/*.phe')
    cancer_names  = list(map(lambda c: icdinfo['cancer'+c] if 'cancer'+c in icdinfo else '',
                         map(lambda path: os.path.splitext(os.path.basename(path))[0], cancer_phenos)))
    make_phe_info(in_phe   = cancer_phenos, 
                  out_path = '/oak/stanford/groups/mrivas/dev-ukbb-tools/phenotypes/extra_info/cancer.tsv.info',
                  name     = cancer_names,
                  field    = 'NA', 
                  table    = 'NA', 
                  basket   = 'NA',
                  app_id   = '16698', 
                  source   = 'NA', 
                  one_file = True)
    # 2. do the same for Rohit's phenotypes
    rh_phenos = ['/oak/stanford/groups/mrivas/private_data/ukbb/16698/phenotypedata/rohit/RH{}.phe'.format(i) for i in range(161)]
    rh_names  = [line.rstrip().split()[1] for line in open('/oak/stanford/groups/mrivas/private_data/ukbb/16698/phenotypedata/rohit/rohit_map.txt', 'r')]
    make_phe_info(in_phe   = rh_phenos, 
                  out_path = '/oak/stanford/groups/mrivas/dev-ukbb-tools/phenotypes/extra_info/rh.tsv.info',
                  name     = rh_names,
                  field    = 'NA', 
                  table    = 'NA', 
                  basket   = 'NA',
                  app_id   = '16698', 
                  source   = 'rohit', 
                  one_file = True)
    # 3. and family history
    fh_phenos = glob.glob('/oak/stanford/groups/mrivas/private_data/ukbb/16698/phenotypedata/familyHistory2/FH*.phe')
    fh_map    = {'FH'+line.split()[0][:4]:line.rstrip().split(None,1)[1].replace(' ','_') for line in open('/oak/stanford/groups/mrivas/private_data/ukbb/16698/phenotypedata/familyHistory2/mapphe.txt', 'r')}
    fh_names  = [fh_map[phe[:6]] if phe[:6] in fh_map else '' for phe in map(os.path.basename,fh_phenos)]
    make_phe_info(in_phe   = fh_phenos,
                  out_path = '/oak/stanford/groups/mrivas/dev-ukbb-tools/phenotypes/extra_info/fh.tsv.info',
                  name     = fh_names,
                  field    = 'NA', 
                  table    = 'NA', 
                  basket   = 'NA',
                  app_id   = '16698', 
                  source   = 'rohit', 
                  one_file = True)
                  

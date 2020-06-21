#!/bin/python
import os
from datetime import date

def make_phe_info(in_phe, out_path, name, field, table, basket, app_id, source, one_file=False):
    # reference 
    pop_dir = '/oak/stanford/groups/mrivas/private_data/ukbb/' + app_id + '/sqc/population_stratification/'
    get_iids = lambda f: set([line.split()[0] for line in open(f, 'r')])
    pops = {pop:get_iids(pop_dir + 'ukb'+ app_id + '_' + pop + '.phe')  
                for pop in ('white_british','non_british_white','african','e_asian','s_asian','semi_related','others')}
    # this is what gets written to file, one item per line below 
    header = "#GBE_ID GBE_NAME FIELD TABLE BASKET APP_ID N N_GBE N_NBW N_AFR N_EAS N_SAS N_SMR N_OTH SOURCE DATE PATH".replace(" ", "\t")
    # in_phe and name are lists of equal sizes -- everything else is constant strings
    if one_file:
        o = open(out_path, 'w')
        o.write(header + '\n')
    # iterate over paths to files and their names (not* GBE IDs)
    for phe_path, gbe_name in zip(in_phe, name):
        # this is a lazy input filtering mechanism from main()
        if not gbe_name: continue
        # load phenotype data 
        if os.path.exists(phe_path):
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
            # get info: each line (approximately) corresponds to one item in the header above
            phe_info = '\t'.join([gbe_id,
                              gbe_name.replace(' ','_'),
                              field,
                              table,
                              basket,
                              app_id,
                              str(sum(map(count,phe.values())))] +
                             [str(sum(map(lambda x:count(x[1]),
                                          filter(lambda x:x[0] in pops[pop], phe.items()))))
                              for pop in ('white_british','non_british_white','african','e_asian','s_asian','semi_related','others')] +                   
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
            # progress bar, effectively: this should be kept in case it's used somewhere else
            print(phe_info[:-1])
            if one_file:
                o.close()
        else:
            print(phe_path + ' not found!')
    return

if __name__ == "__main__":
    #return
    #import glob
    import pandas as pd
    phenotype_info = pd.read_csv('../../05_gbe/phenotype_info.tsv', sep='\t', dtype=object, keep_default_na=False)
    all_phenos = list(phenotype_info['PATH'])
    all_names = list(phenotype_info['GBE_NAME'])
    all_fields = list(phenotype_info['FIELD'])
    all_tables = list(phenotype_info['TABLE'])
    all_baskets = list(phenotype_info['BASKET'])
    all_sources = list(phenotype_info['SOURCE'])
    for phe, gbe_name, field_id, table_id, basket_id, source_id in zip(all_phenos, all_names, all_fields, all_tables, all_baskets, all_sources):
        make_phe_info(in_phe   = [phe],
                      out_path = os.path.join(os.path.dirname(os.path.dirname(phe)), "info"),
                      name     = [gbe_name],
                      field    = field_id, 
                      table    = table_id,
                      basket   = basket_id,
                      app_id   = '24983',
                      source   = source_id, 
                      one_file = False)

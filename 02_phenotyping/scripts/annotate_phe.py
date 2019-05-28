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
        # progress bar, effectively: this should be kept in case it's used somewhere else
        print(phe_info[:-1])
    if one_file:
        o.close()
    return


if __name__ == "__main__":
    import glob
    import pandas as pd
    all_phenos = []
    all_names = []
    all_fields = []
    all_tables = []
    all_baskets = []
    all_sources = []
    meta_table = pd.read_table('../tables/gbe_sh_input_params.tsv')
    map_tables = []
    for gbe_id_col, field_id_col, desc_col, table_path in zip(meta_table['nameCol (GBE ID)'], meta_table['fieldCol (field_ID)'], meta_table['descCol'], meta_table['Table path']):
        colnames = [x[1] for x in sorted(zip([gbe_id_col, field_id_col, desc_col], ['GBEID', 'FieldID', 'Name']), key=lambda i:i[0])]
        map_table = pd.read_table('../tables/' + table_path, usecols=[gbe_id_col, field_id_col, desc_col], sep='\t', skiprows=1, dtype=str, names=colnames)
        map_table['Source'] = table_path
        map_tables.append(map_table)
    complete_table = pd.concat(map_tables)
    id_to_field_name = dict(zip(complete_table['GBEID'], zip(complete_table['Name'], complete_table['FieldID'], complete_table['Source'])))
    for basket in ['2001440']:
        phenos = glob.glob('/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/' + basket + '/current/phe/*.phe')
        phenos = [os.path.realpath(pheno) for pheno in phenos]
        gbe_ids = [os.path.splitext(os.path.basename(path))[0] for path in phenos]
        tables = ['9797' if gbe_id.startswith('MED') else os.path.basename(os.path.normpath(os.path.dirname(os.path.dirname(pheno)))) for gbe_id, pheno in zip(gbe_ids, phenos)]
        with open('../../05_gbe/icdinfo.txt', 'r') as new_icd:
            new_icdinfo = {line.split()[0]:line.split()[2] for line in new_icd}
        names = [id_to_field_name[gbe_id][0] if gbe_id in id_to_field_name else new_icdinfo[gbe_id] for gbe_id in gbe_ids]
        fields = ['20003' if gbe_id.startswith('MED') else id_to_field_name[gbe_id][1] for gbe_id in gbe_ids]
        sources = ['NA' if gbe_id.startswith('MED') else id_to_field_name[gbe_id][2] for gbe_id in gbe_ids]
        tables = ['9797' if gbe_id.startswith('MED') else os.path.basename(os.path.normpath(os.path.dirname(os.path.dirname(pheno)))) for pheno in phenos]
        baskets = [basket for gbe_id in gbe_ids]
        all_phenos.extend(phenos)
        all_names.extend(names)
        all_fields.extend(fields)
        all_tables.extend(tables)
        all_baskets.extend(baskets)
        all_sources.extend(sources)

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

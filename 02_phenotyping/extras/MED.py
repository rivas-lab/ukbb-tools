
import numpy as np
import pandas as pd
import glob, os

# program these
data_root = '/oak/stanford/groups/mrivas/private_data/ukbb/24983/'
out_dir   = '/oak/stanford/groups/mrivas/dev-ukbb-tools/phenotypes/9796/20003/'

# get input data
src = glob.glob(os.path.join(data_root, 'phenotypedata/download/9796/newest/ukb*.tab'))[0]
input_data = pd.read_table(src, usecols=lambda x:'20003' in x.split('.') or x == 'f.eid', index_col='f.eid', dtype='str')
phe_ref = dict([line.split() for line in open('MED_input.tsv')])

# useful for phenotype info
header="#GBE_ID GBE_NAME FIELD TABLE BASKET APP_ID N N_GBE N_AFR N_EAS N_SAS SOURCE DATE PATH".replace(" ", "\t")
popref = {pop:pd.read_table(os.path.join(data_root, 'sqc/population_stratification/ukb24983_{}.phe'.format(pop)),
                            index_col=0).iloc[:,0].to_dict() for pop in ('white_british', 'african', 'e_asian', 's_asian')}
          

# iterate through reference info
for i,(code,name) in enumerate(phe_ref.items()):
    # ok what are we doing
    new_phe = input_data.apply(lambda x:code in x.tolist(), axis=1).to_dict()
    out_phe = os.path.join(out_dir, 'MED'+code+'.phe')
    # write phenotype file
    with open(out_phe, 'w') as o:
        o.write("\n".join(map(lambda x:' '.join((str(x[0]),str(x[0]),'2' if x[1] else '1')), new_phe.items())))
    # write phenotype info
    with open(out_phe+'.info', 'w') as o:
        o.write(header+'\n')
        o.write('\t'.join(['MED'+code, 
                           name, 
                           '20003', 
                           os.path.splitext(os.path.basename(src))[0].replace('ukb',''), 
                           '9796', 
                           '24983', 
                           str(sum(new_phe.values()))] + 
                          [str(len(set(popref[p].keys()) & set([k for k,v in new_phe.items() if v]))) for p in 
                           ['white_british', 'african', 'e_asian', 's_asian']] + 
                          ['MED.py', 
                           str(pd.Timestamp.today()).split()[0], 
                           out_phe]
                         ) + '\n') 

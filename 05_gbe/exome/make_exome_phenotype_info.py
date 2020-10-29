import os
import pandas as pd
import glob

# all exome
print("Reading in exome psam file...")
exome = {}
exome['all'] = pd.read_csv("/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/ukb24983_exomeOQFE.psam", sep='\t', dtype=str)

# get the rest of the populations
print('Reading in population files...')
for pop in ['african', 'e_asian', 'non_british_white', 's_asian', 'white_british', 'semi_related', 'others']:
    exome[pop] = exome['all'].merge(pd.read_csv("/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_" + pop + ".phe", sep='\s+', header=None, usecols=[1], names=['IID'], dtype=str), on='IID')

# do the thing
print("Progressing through phenotypes...")
with open("../array-combined/phenotype_info.tsv", "r") as f, open("200k/exome_phenotype_info.tsv","w") as o:
    for n,line in enumerate(f):
        # write header when we encounter the header
        if n == 0: 
	    o.write(line)
            continue
        # ok go
        phe_info = line.rstrip().split('\t')
        phe_data = pd.read_csv(phe_info[-1], sep='\s+', header=None, names=['#FID','IID','PHENO'], dtype='str')
        # if all values are 1,2,-9 then the trait is binary, otherwise QT
        if all([i in ['-9','1','2'] for i in phe_data['PHENO'].value_counts().index.tolist()]):
            query = 'PHENO == "2"'
        else:
            query = 'PHENO != "-9"'
        # this too shall pass
        pop_count = {}
        for pop,pop_data in exome.items():
            # count up individuals
            pop_count[pop] = str(phe_data.query(query).merge(pop_data, on='IID').shape[0])
        # write to file 
        o.write("\t".join(phe_info[:6] + [pop_count[p] for p in ['all','white_british','non_british_white','african','e_asian','s_asian','semi_related','others']] + phe_info[-3:]) + "\n")

epi = pd.read_table('200k/exome_phenotype_info.tsv')
short = pd.read_table('../array-combined/icdinfo.shortnames.tsv')

merged = short.merge(epi, how='left', left_on='GBE_ID', right_on='#GBE_ID')

merged = merged[['GBE_category', 'GBE_ID', 'N_GBE', 'GBE_NAME_x', 'GBE_short_name', 'GBE_short_name_len', 'Units_of_measurement']]
merged.columns = ['GBE_category', 'GBE_ID', 'GBE_N_EXOME', 'GBE_NAME', 'GBE_short_name', 'GBE_short_name_len', 'Units_of_measurement']
merged.to_csv('200k/icdinfo.shortnames.exome.tsv', index=False, sep='\t')        

import os
import pandas as pd
import glob

# all exome
exome = {}
exome['all'] = pd.read_table("/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/spb/data/ukb_exm_p_spb.psam", dtype=str)

# get the rest of the populations
for f in glob.glob("/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_*.phe"):
    pop = os.path.basename(f).replace('.phe','').replace('ukb24983_','')
    exome[pop] = exome['all'].merge(pd.read_table(f, header=None, usecols=[1], names=['IID'], dtype=str), on='IID')

# do the thing
with open("phenotype_info.tsv", "r") as f, open("exome_phenotype_info.tsv","w") as o:
    for n,line in enumerate(f):
        # write header when we encounter the header
        if n == 0: 
	    o.write(line)
            continue
        # ok go
        phe_info = line.rstrip().split('\t')
        phe_data = pd.read_table(phe_info[-1], header=None, names=['#FID','IID','PHENO'], dtype='str')
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
        o.write("\t".join(phe_info[:6] + [pop_count[p] for p in ['all','white_british','african','e_asian','s_asian']] + phe_info[-3:]) + "\n")
        

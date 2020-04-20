#!/bin/python
from __future__ import print_function
from datetime import date

_README_='''
README:

A script to create a "master" phenotype file from resources in the ukbb-tools repo.

Author: Matthew Aguirre (SUNET: magu)
'''

# what do?
out_file = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.{}.phe'.format(str(date.today()).replace('-',''))

# get list of phenotypes
phe_in = {}
with open('../05_gbe/phenotype_info.tsv', 'r') as f, open(out_file+'.info.tsv', 'w') as o:
    header = "#GBE_ID GBE_NAME FIELD TABLE BASKET APP_ID N N_NBW N_GBE N_AFR N_EAS N_SAS SOURCE DATE PATH".replace(" ", "\t")
    o.write(header + '\n')
    for i,line in enumerate(f):
        info = line.rstrip().split('\t')
        # process header
        if i == 0: 
            cols = {field:col for col,field in enumerate(info)}
        # load data
        else:
            if int(info[cols['N_GBE']]) >= 100: 
                phe_in[info[cols['#GBE_ID']]] = {n:info[cols[n]] for n in ['APP_ID','PATH']}
                o.write(line)
            
# load up individual data map: all IDs in the output will be for app 24983
to_24983 = {}
with open('/oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_16698_mapping.tsv', 'r') as f:
    for line in f:
        l = line.split()
        to_24983[l[1]] = l[0]

# get data
inds    = {} # phenotype data, keyed by iid
phe_out = [] # names

# start with the gwas covariates
with open('/oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar.phe') as f:
    for i,line in enumerate(f):
        l = line.rstrip().split('\t')
        if i == 0:
            phe_out += l[2:]
        else:
            inds[l[0]] = l[2:]
n_covar = len(phe_out)

print("Loading {0} phenotypes and {1} covariates:".format(len(phe_in), n_covar))

# now do the phenotype data
for phe,app,path in map(lambda (k,v): (k,v['APP_ID'],v['PATH']),phe_in.items()):
    # keep track of phenotype, then import
    phe_out.append(phe)
    # don't assume the individuals are the same across files
    inds_not_seen = set(inds.keys()) 
    # the above will be pruned as we walk though the file
    with open(path, 'r') as f:
        for line in f:
            ind,_,value = line.split()
            if app == '16698' and ind in to_24983:
                ind = to_24983[ind]
            if ind not in inds:
                inds[ind] = ['-9' for _ in range(len(phe_out)-1)] + [value]
            else:
                inds[ind].append(value)
                inds_not_seen.remove(ind)
        # add empty values for missing individuals
        for ind in inds_not_seen:
            inds[ind].append('-9')
    print("In: {1}/{2} ({0})".format(phe, len(phe_out), len(phe_in)+n_covar))

# write to file
sep='\t'
print("writing to file...")
with open(out_file, 'w') as f:
    f.write(sep.join(['FID', 'IID'] + phe_out) + '\n')
    for ind in sorted(inds.keys()):
        f.write(sep.join([ind, ind] + inds[ind]) + '\n')

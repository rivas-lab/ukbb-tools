#!/bin/python
import glob
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.decomposition import TruncatedSVD

# pvalue cutoff for var::phe association
p = 0.001
# z-statistic instead of beta?
z = False
# z-score values by phenotype?
c = True
# components for tsvd
n = 100

# find files to import
files = glob.glob('../summary_stats/*.gz')

# only keep variants not in LD, MAD > 0.01%, and QC'd
qc_file = '/oak/stanford/groups/mrivas/private_data/ukbb/variant_filtering/variant_filter_table.tsv.gz'
qc      = pd.read_table(qc_file, usecols=['ID','ALT','maf','ld_indep','all_filters'], index_col='ID')
var_set = set(qc[(qc['maf'] > 0.0001) & (qc['all_filters'] == 0) & (qc['ld_indep'] == True)].index.tolist())
var_alt = qc.loc[list(var_set),'ALT'].tolist()

# only keep traits with n>100 cases or observations
with open('../summary_stats/train_phe_counts.tsv', 'r') as f:
	phe_n = {line.split()[0]:int(line.rstrip().split()[1]) for line in f}

# load reference for rsid to minor allele

# helper function to load files
def get_summary_stats(f, p=0.01, z=True):
	# load data, perform p-value and LD pruning
	maxerr = 0.8 if 'linear' in f else 0.2 if 'logistic' in f else 0
	x = pd.read_table(f, index_col=2) 
	x = x[(x['P'] < p) & (x['TEST'] == "ADD") & (x.index.isin(var_set)) & (x['SE'] < maxerr)]
	# take LOR if binary, manage z-scoring
	if 'linear' in f:
		return x['BETA'].divide(x['SE'] if z else 1)
	elif 'logistic' in f:
		return x['OR'].apply(np.log).divide(x['SE'] if z else 1)

# load data
data = pd.SparseDataFrame(index=var_set)
for f in files:
	phe_id = f.split('/')[-1].split('.')[1]
	if phe_n[phe_id] > 100:
		data[phe_id] = get_summary_stats(f,p,z)
if c:
	data = data.subtract(data.mean()).divide(data.std())

# some phenotypes are duplicates
def dup_check(x, phe_list=data.columns):
	if len(x) <= 5:
		return False
	if 'INI100' in x:
		return ('INI'+x[-4:] in phe_list) or ('INI'+x[-5:] in phe_list)
	if any(['HC{}00'.format(i) in x for i in range(1,19)]):
		return ('BIN'+x[-4:] in phe_list) or ('BIN'+x[-5:] in phe_list)

data = data.drop(labels=filter(dup_check, data.columns), axis=1)

# name and save
dataset_name = '_'.join(('all',
                         'z' if z else 'beta',
                         'center' if c else 'nonCenter',
                         'p'+str(p).replace('.',''),
                         str(n)+'PCs',
                         str(pd.Timestamp.today()).split()[0].replace('-','')))

data.to_pickle(dataset_name + '.full_df.pkl.gz', compression='gzip')

# do the analysis
matt = TruncatedSVD(n_components=n, n_iter=20, random_state=24983)
US = matt.fit_transform(csr_matrix(data.fillna(value=0).values))

# save the results
np.savez(dataset_name,
         U = US/matt.singular_values_,
         V = matt.components_.T,
         D = matt.singular_values_,
         variance_explained = matt.explained_variance_,
         variance_explained_ratio = matt.explained_variance_ratio_,
         label_phe_code = np.array(data.columns),
         label_var = np.array(data.index),
         label_var_minor_allele = np.array(var_alt) 
) 


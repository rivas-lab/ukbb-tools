import argparse
import glob
import collections
import gzip
import os
import pdb
import subprocess
import sys

import numpy as np
import pandas as pd

def read_plink2_sumstats(fn, regression_type, header=True):
    """Read a plink 2 output file into a pandas DataFrame.

    Parameters
    ----------
    fn : str 
        Path to the plink file. The file can be gzipped or not.

    regression_type : str
        This needs to be either 'logistic' or 'linear'.

    header : str
        True if the file has a header (this is generally the case unless the
        file has been processed after it was created).

    Returns
    -------
    res : pandas.DataFrame
        Dataframe with results.
    """

    dtypes_dict = dict()
    dtypes_dict['logistic'] = collections.OrderedDict([
	('#CHROM', str), ('POS', int), ('ID', str), ('REF', str), ('ALT1', str),
        ('FIRTH?', str), ('TEST', str), ('OBS_CT', int), ('OR', float),
	('SE', float), ('T_STAT', float), ('P', float)
    ])
    dtypes_dict['linear'] = collections.OrderedDict([
	('#CHROM', str), ('POS', int), ('ID', str), ('REF', str), ('ALT1', str),
	('TEST', str), ('OBS_CT', int), ('BETA', float),
	('SE', float), ('T_STAT', float), ('P', float)
    ])

    assert(regression_type in set(dtypes_dict.keys()))
    dtypes = dtypes_dict[regression_type]

    if header is None:
        if fn[-3:] == '.gz':
            with gzip.open(fn, 'r') as f:
                line = f.readline()
	else:
            with open(fn, 'r') as f:
                line = f.readline()
        header = line [0] == '#'
    try:
        if header:
            res = pd.read_table(fn, index_col=2, dtype=dtypes, low_memory=False)
        else:
	    cols = list(dtypes.keys())
            res = pd.read_table(fn, index_col=2, dtype=dtypes, names=cols,
                                low_memory=False)
    except pd.io.common.EmptyDataError:
        sys.stderr.write('No data in {}\n'.format(fn))
        sys.exit()

    res.columns = [x.replace('#', '') for x in res.columns]
    res.columns = [x.replace('ALT1', 'ALT') for x in res.columns]
    return(res)

def read_logistic2(fn, header=True):
    """Read a plink 2 output file of type glm.logistic into a pandas DataFrame.

    Parameters
    ----------
    fn : str 
        Path to the plink file. The file can be gzipped or not.

    header : str
        True if the file has a header (this is generally the case unless the
        file has been processed after it was created).

    Returns
    -------
    res : pandas.DataFrame
        Dataframe with results.
    """
    return(read_plink2_sumstats(fn, regression_type='logistic', header=header))

def read_linear2(fn, header=True):
    """Read a plink 2 output file of type glm.linear into a pandas DataFrame.

    Parameters
    ----------
    fn : str 
        Path to the plink file. The file can be gzipped or not.

    header : str
        True if the file has a header (this is generally the case unless the
        file has been processed after it was created). False if no header. Pass
        None if it's unknown whether the file has a header.

    Returns
    -------
    res : pandas.DataFrame
        Dataframe with results.

    """
    return(read_plink2_sumstats(fn, regression_type='linear', header=header))

def filter_linear_output(res):
    """Filter the PLINK 2 linear output"""
    res = res[res['TEST'] == 'ADD']
    res = res.dropna(subset=['BETA'])
    res['SE'] = res['SE'] 
    # dfs.append(res[['BETA', 'SE', 'P']])
    return(res)

def filter_logistic_output(res):
    """Filter the PLINK 2 logistic output"""
    res = res[res['TEST'] == 'ADD']
    res = res.dropna(subset=['OR'])
    res['BETA'] = np.log(res['OR'])
    res['SE'] = res['SE'] 
    # dfs.append(res[['BETA', 'SE', 'P']])
    return(res)
            
def get_linear_output(fns, header=False):
    """Read plink linear output files, combine them, and
    filter/transform as needed."""
    dfs = []
    for fn in fns:
        res = read_linear2(fn, header=None)
	dfs.append(filter_linear_output(res))
    return pd.concat(dfs)

def get_logistic_output(fns, header=False):
    """Read plink logistic output files, combine them, and
    filter/transform as needed."""
    dfs = []
    for fn in fns:
        res = read_logistic2(fn, header=None)
	dfs.append(filter_logistic_output(res))
    return(pd.concat(dfs))

def get_results(code):
    if code[0:2] == 'HR':
        files = [('/oak/stanford/groups/mrivas/private_data/ukbb/16698/cal/gwas/highconfidenceqc_nhs_newest/'
                  'ukb16698_v2.HC{}.PHENO1_c{}.glm.logistic.hybrid.gz'.format(
                      code[2:], x)) for x in range(1, 23)]
        code_data = get_logistic_output(files)
    elif code[0:2] == 'VQ':
        files = [('/oak/stanford/groups/mrivas/private_data/ukbb/16698/cal/gwas/highconfidenceqc_verbal_newest/'
                  'ukb16698_v2.HC{}.PHENO1_c{}.glm.logistic.hybrid.gz'.format(
                      code[2:], x)) for x in range(1, 23)]
        code_data = get_logistic_output(files)
    else:
        dy = '/oak/stanford/groups/mrivas/private_data/ukbb/gwas_results'
        gwas_info = pd.read_table(os.path.join(dy, 'hybrid.tsv'))
        tdf = gwas_info[gwas_info.code == code]
        # For now, I'm only supporting phenotypes with results split by chromosome.
        assert tdf.shape[0] == 22, 'Wrong amount of files for {}'.format(code)
        files = [os.path.join(dy, 'hybrid', x) for x in tdf['fn']]
        if tdf['regtype'].values[0] == 'logistic':
            code_data = get_logistic_output(files)
        elif tdf['regtype'].values[0] == 'linear':
            code_data = get_linear_output(files)
    return(code_data)

def _read_ld_scores(path):
    fns = glob.glob(os.path.join(path, '*ldscore.gz'))
    lds = [pd.read_table(fn) for fn in fns]
    lds = pd.concat(lds)
    lds.index = lds.CHR.astype(str) + ':' + lds.BP.astype(str)
    se = pd.Series(lds.index)
    # Some (chromosome, position) pairs are not unique on the UKB array. This
    # causes some headaches, so I'm just going to remove those ones here (which
    # will cause them to be removed from the summary stats too).
    vc = se.value_counts()
    remove = vc[vc > 1].index
    lds = lds.drop(remove)
    return(lds)

def _merge_with_ld_scores(res, lds):
    # Make sure summary stats and LD scores have same SNP names.
    res.index = res['CHROM'].astype(str) + ':' + res.POS.astype(str)
    shared = set(res.index) & set(lds.index)
    res = res.loc[shared]
    res.loc[res.index, 'ID'] = lds.loc[res.index, 'SNP'].values
    return(res)

def make_file(pheno, ld_scores, outdir, keep, fn_logistic = None, fn_linear = None):
    assert(fn_logistic is None or fn_linear is None)
    if(fn_logistic is not None):
	files = glob.glob(fn_logistic)
        data = get_logistic_output(files)
    elif(fn_linear is not None):
	files = glob.glob(fn_linear)
        data = get_linear_output(files)
    else:
        data = get_results(pheno)
    
    # Filter out variants that we don't have LD scores for.
    lds = _read_ld_scores(ld_scores)
    data = _merge_with_ld_scores(data, lds)
    # Filter further if needed.
    if keep:
        keep = set(pd.read_table(keep, header=None, squeeze=True))
        data.index = data['ID']
        data = data.loc[set(keep) & set(data.index)]
    # Write output file.
    data.to_csv(os.path.join(outdir, '{}_ldsc.tsv'.format(pheno)),
                sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(
        description=('Process plink output for use with LD score regression.'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('phenotype_code', 
                        help='Phenotype code for the trait of interest.')
    parser.add_argument('-o', metavar='outdir', default='.',
                        help=('Output directory for files.'))
    ld_scores = os.path.realpath(os.path.join(os.path.split(__file__)[0],
                                              '../private_output/ukb_ld_scores/TWINSUK'))
    parser.add_argument('-ld', metavar='ld_scores', default=ld_scores,
                        help=('Directory with LD scores.'))
    keep = os.path.realpath(os.path.join(os.path.split(__file__)[0],
                                         '../misc/variants_to_use_ldsc.txt'))
    parser.add_argument('-k', metavar='keep', default=keep,
                        help=('Only use variants present in this file for '
                              'analysis.'))
    file_help=(
	'Input file of the PLINK summary statistics ({}).'
	'You may use glob. If you use this option, '
	'we will not check the number of files (chromosomes).'
    )
    parser.add_argument('--logistic_file', metavar='plink2_sumstats_file', 
			help=file_help.format('logistic'))
    parser.add_argument('--linear_file', metavar='plink2_sumstats_file', 
			help=file_help.format('linear'))
    args = parser.parse_args()

    pheno = args.phenotype_code
    ld_scores = os.path.realpath(args.ld)
    outdir = args.o
    keep = os.path.realpath(args.k)
    fn_logistic = args.logistic_file
    fn_linear   = args.linear_file
    assert(fn_logistic is None or fn_linear is None)
    make_file(pheno, ld_scores, outdir, keep, fn_logistic, fn_linear)

if __name__ == '__main__':
    main()

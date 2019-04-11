from __future__ import print_function

_README_ = '''
-------------------------------------------------------------------------
phewas.py

Query SciDB (GBE) and get the PheWAS sumstats

Author: Yosuke Tanigawa (ytanigaw@stanford.edu)
Date: 2019/04/11
-------------------------------------------------------------------------
'''

import os, sys, collections, argparse
from scidbpy import connect
import numpy as np
import pandas as pd

# import some relevant functions from GBE repo
gbe_repo = os.path.join('/', 'opt', 'biobankengine', 'GlobalBioBankEngineRepo')
sys.path.append(os.path.join(gbe_repo, 'gbe_browser'))
import config, lookups


def variant_page_data_prep_sub(icdstats, sort_key='log10pvalue'):
    plot_d_raw = collections.defaultdict(list)
    keys = icdstats[0].keys()
    for key in keys:
        plot_d_raw[key] = np.array([x[key] for x in icdstats])
    plot_df = pd.DataFrame(plot_d_raw)
    #.sort_values(
    #    by=['Group', sort_key], ascending=[True, False]
    #)
    plot_d_dict = collections.defaultdict(collections.defaultdict)
    
    groups = sorted(set(plot_df['Group']))
    for group in groups:
        for key in keys:
            plot_d_dict[group][key] = list(plot_df[plot_df['Group'] == group][key])
    for group in groups:
        for key in ['OR', 'LOR', 'L95OR', 'U95OR', 'pvalue', 'SE', 'log10pvalue']:
            plot_d_dict[group][key] = [float(x) for x in plot_d_dict[group][key]]
    for group in groups:
        #error_bar = {'L95OR': -1, 'U95OR': 1}
        #for key in error_bar.keys():
        #    diff = np.array(plot_d_dict[group][key]) - np.array(plot_d_dict[group]['LOR'])
        #    plot_d_dict[group]['d{}'.format(key)] = [0 if np.isnan(x) else np.abs(x) for x in diff]
        plot_d_dict[group]['196SE'] = list( 1.96 * np.array(plot_d_dict[group]['SE']) )

    for group in groups:
        if group in set(['INI', 'INI_FC', 'BROADQT']):
            beta_or_lor = 'BETA'
            beta_or_lor_val = plot_d_dict[group]['LOR']
            beta_or_lor_l95 = np.array(plot_d_dict[group]['LOR']) - 1.96 * np.array(plot_d_dict[group]['SE'])
            beta_or_lor_u95 = np.array(plot_d_dict[group]['LOR']) + 1.96 * np.array(plot_d_dict[group]['SE'])

        else:
            beta_or_lor = 'OR'
            beta_or_lor_val = np.exp(np.array(plot_d_dict[group]['LOR']))
            beta_or_lor_l95 = np.exp(np.array(plot_d_dict[group]['LOR']) - 1.96 * np.array(plot_d_dict[group]['SE']))
            beta_or_lor_u95 = np.exp(np.array(plot_d_dict[group]['LOR']) + 1.96 * np.array(plot_d_dict[group]['SE']))

        group_len = len(plot_d_dict[group]['icd'])
        plot_d_dict[group]['text'] = [
            '{}. Case: {}, P-value: {:.3e}, {} = {:.5f} (95% [{:.5f}, {:.5f}]), SE = {:.5f}'.format(
                ''.join([c if c != '_' else ' ' for c in x[0]]), x[1], x[2], x[3], x[4], x[5], x[6], x[7]
            ) for x in zip(
                plot_d_dict[group]['Name'],
                plot_d_dict[group]['Case'],
                plot_d_dict[group]['pvalue'],
                [beta_or_lor] * group_len,
                beta_or_lor_val,
                beta_or_lor_l95,
                beta_or_lor_u95,
                plot_d_dict[group]['SE'],
                #plot_d_dict[group]['L95OR'],
                #plot_d_dict[group]['U95OR'],
            )
        ]
    return plot_d_dict

def get_data(variant_str, db):

    chrom, pos = variant_str.split('-')
    icdstats = lookups.get_icd_by_chrom_pos(db, chrom, pos)

    indexes = []
    seend = {}

    for idx in range(0,len(icdstats)):
        # ICD10=T81/Name=Complications of procedures, not elsewhere classified/Chapter=T/L95OR=0.97/U95OR=2.04/OR=1.40/pvalue=0.0756463/l10pval=1.12/Case=1229
        item = icdstats[idx]
        icd10 = item['icd']
        item['Code'] = icd10
        # icd10info = lookups.get_icd_info(db, icd10)
        if 'Name' not in item:
            item['Name'] = 'NA'
            item['Group'] = 'NA'
            item['OR'] = 1
            item['LOR'] = 0
            item['L95OR'] = 1
            item['U95OR'] = 1
            item['pvalue'] = 1
            item['l10pval'] = 0
            item['Case'] = 'NA'
            item['SE'] = 0
            indexes.append(idx)
        else:
            # item['Name'] = icd10info[0]['Name']
            item['Group'] = icd10[0] # default value
            groups = ['RH', 'FH', 'HC', 'cancer', 'ADD', 'INI', 'MED', 'BIN', 'BRMRI', 'BROADBIN', 'BROADQT', 'INI_FC', 'BIN_FC']
            for group in groups:
                if icd10.startswith(group):
                    item['Group'] = group
                    break
            item['OR'] = format(float(item['or_val']), '.4g')
            item['LOR'] = format(float(item['lor']), '.4g')
            item['L95OR'] = format(float(item['l95or']), '.4g')
            item['U95OR'] = format(float(item['u95or']), '.4g')
            item['pvalue'] = format(float(item['pvalue']), '.4g')
            item['l10pval'] = format(float(item['log10pvalue']), '.4g')
            item['SE'] = format(float(item['se']), '.4g')
            if float(item['pvalue']) == 0:
                item['pvalue'] = numpy.finfo(float).eps
                item['pvalue'] = format(float(item['pvalue']),'.4g')
                item['l10pval'] = 250
            # item['Case'] = icd10info[0]['Case']
            se =  format(float(item['se']), '.4g')
            if (
		(float(item['l10pval']) < 1) or 
		(float(se) >= .5) or 
		(float(se) >= .08 and item['OR'] == item['LOR']) or 
		(int(item['Case']) <= 100)  or 
		(item['Code'] == "HC67") or 
		(icd10 in seend)
	    ):
                indexes.append(idx)
            seend[icd10] = icd10


    d = variant_page_data_prep_sub(icdstats)
    return d

def filter_PheWAS_by_p_value(data_dict, p_val_thr):     
    new_d = collections.OrderedDict()
    for category, cat_data in data_dict.items():
        data_filter = np.array(
            [x <= p_val_thr for x in cat_data['pvalue']]
        )
        new_d[category] = collections.OrderedDict([
            (k, np.array(v)[data_filter]) 
            for k, v in zip(cat_data.keys(), cat_data.values())
        ])
    return new_d        

def df_PheWAS(data_dict):    
    if(len(data_dict) == 0):
        return pd.DataFrame()
    else:
        float_cols = [
            'l10pval', 'log10pvalue',
            '196SE',
            'lor', 'LOR', 'or_val',
            'l95or', 'u95or',             
            'L95OR', 'U95OR',            
        ]
        df_concat = pd.DataFrame()
        for data_per_category in data_dict.values():
            df_concat = pd.concat(
                [
                    df_concat, 
                    pd.DataFrame(data_per_category)
                ], 
                ignore_index=True
            )
        for col in float_cols:
            df_concat[col] = np.array([
                float(x) for x in df_concat[col]
            ])
        df_concat['is_binary'] = df_concat['Group'].map(
            lambda x: x not in set(['INI', 'INI_FC', 'BROADQT'])
        )            
        return df_concat

def df_PheWAS_filter_by_case_count(df, min_case_count):
    return df[df['Case'].map(lambda x: x >= min_case_count)]

def get_data_as_data_frame(query, db, p_value_thr = 0.001, min_case_count = 100, phenotype_list=None):
    query_res = get_data(query, db)
    res_tbl = df_PheWAS_filter_by_case_count(
        df_PheWAS(filter_PheWAS_by_p_value(query_res, p_value_thr)),
        min_case_count
    )[['affyid', 'Code', 'l10pval', 'LOR', 'SE', 'chrom', 'pos', 'Name']]
    if phenotype_list is None:
	return(res_tbl)
    else:
	return(res_tbl[ res_tbl['Code'].map(lambda x: x in set(phenotype_list)) ])

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=_README_
    )
    parser.add_argument('-i', metavar='i', required=True,
                        help='query (comma-separated coordinates like 5-145895394,11-14865399)')    
    parser.add_argument('-o', metavar='o', required=True,
                        help='output npz file')
    parser.add_argument('--port', metavar='P', default=8080,
                        help='SciDB port (default: 8080)')
    parser.add_argument('-n', metavar='N', default=100,
                        help='minimum case count (default: 100)')
    parser.add_argument('-p', metavar='p', default=0.001,
                        help='p-value threshold (default: 0.001)')
    parser.add_argument('-l', metavar='L', default=None,
                        help='a file that has a list of phenotypes (GBE_IDs)')
  
    args = parser.parse_args()
    db = connect('http://localhost:{}'.format(args.port))

    if(args.l is not None):
	with open(args.l) as f:
	    phe_list = [x for x in f.read().splitlines() if len(x) > 0]
    else:
	phe_list = None

    query_res = pd.concat([ 
	get_data_as_data_frame(
	    i, db, 
	    p_value_thr = float(args.p), 
	    min_case_count = int(args.n), 
	    phenotype_list= phe_list
	) for i in list(set(args.i.split(','))) 
    ])
    query_res.to_csv(args.o, sep='\t', index=False)

if __name__ == '__main__':
    main()


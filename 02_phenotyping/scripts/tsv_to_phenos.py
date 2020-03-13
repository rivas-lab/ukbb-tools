#!/bin/python
import os,glob
import pandas as pd
from make_phe import *
from annotate_phe import make_phe_info

_README_="""
This script is designed to process table-defined phenotypes (from a Rivas lab computing session) into *.phe files usable for downstream analysis with PLINK and similar tools. 

To define all phenotypes from a table, identify the (zero-indexed!) columns specified below, and pass them to the script. If you would like to define only one phenotype from the table, you can do that by adding the flag --only-this-row, followed by the row (zero-indexed, such that the first non-header row is number 0) you'd like to run.

Author: Matthew Aguirre and Guhan Venkataraman (SUNETs: magu and guhan)
Updated: 2020/01/29
"""

def get_phe_definitions(in_tsv, field_col, name_col, case_col, ctrl_col, excl_col, qtfc_col, 
             desc_col, ann_col, ann_d_col, ftbd, header=True, make_this_one=None):

    # get phenotype definitions from an input table, with option to only extract one row
    phe_info = {}
    if make_this_one is not None:
        if isinstance(make_this_one, list):
            make_this_one = int(make_this_one[0])
        else:
            make_this_one = int(make_this_one)
    with open(in_tsv, 'r') as f:
        for i, line in enumerate(f):
            if header and i==0:
                if make_this_one is not None:
                    make_this_one += 1
                continue
            if make_this_one is not None and i != make_this_one:
                continue
            fields = line.rstrip().split('\t')
            if not ((fields[name_col] and fields[desc_col] and fields[field_col]) and (fields[case_col] or fields[ctrl_col] or fields[excl_col] or fields[qtfc_col] or fields[ann_col] or fields[ann_d_col])): 
                if fields[name_col]:
                    print("No annotation data present for " + fields[name_col] + ".")
                continue
            field_metadata = ftbd[ftbd['Field_ID'] == fields[field_col]]
            # the if switch is in case empty cells get pruned from the TSV
            phe_info[ fields[name_col] ] = {'case':      fields[case_col] if case_col < len(fields) else '',
                                            'control':   fields[ctrl_col] if ctrl_col < len(fields) else '',
                                            'qt_order':  fields[qtfc_col] if qtfc_col < len(fields) else '',
                                            'exclude':   fields[excl_col] if excl_col < len(fields) else '',
                                            'table_id':  field_metadata['Table_ID'].values[0] if len(field_metadata) > 0 else '',
                                            'basket_id': field_metadata['Basket_ID'].values[0] if len(field_metadata) > 0 else '',
                                            'field_id':  fields[field_col] if field_col < len(fields) else '',
                                            'desc':      fields[desc_col] if desc_col < len(fields) else ''}
    return phe_info


def define_phenos(in_tsv, field_col, name_col, case_col, ctrl_col, 
                  excl_col, qtfc_col, desc_col, ann_col, ann_d_col, ftbd,
                  header=True, all_ctrl=False, make_this_one=None):
    home_dir='/oak/stanford/groups/mrivas/private_data/ukbb/24983/phenotypedata'
    # iterate over phenotypes in the input file
    for phe_name, phe_values in get_phe_definitions(in_tsv, field_col, name_col, case_col, ctrl_col, excl_col, 
                qtfc_col, desc_col, ann_col, ann_d_col, ftbd, header, make_this_one).items():
        phe = os.path.join(home_dir, phe_values['basket_id'], phe_values['table_id'], 'phe', phe_name + '.phe')
        log = os.path.join(os.path.dirname(os.path.dirname(phe)), "logs/{0}.log".format(phe_name))
        info_dir = os.path.join(os.path.dirname(os.path.dirname(phe)), "info")
        info_file = os.path.join(info_dir, "{}.info".format(phe_name))
        #if os.path.exists(phe) and os.path.exists(info_file):
        #    print("Phenotype for " + phe_name + " already exists at " + phe + ".")
        #else:
        #    print(phe_values)
        if True:
            tab_f = 'ukb{}.tab'.format(phe_values['table_id'] )
            # this will throw an indexing error if a bad table is supplied
            tsv = os.path.join(home_dir, phe_values['basket_id'], phe_values['table_id'], 'download', tab_f)
            # assume binary if we have a case definition, else assume qt
            if phe_values['case']: 
                # this and create_qt_phe_file below are implemented in make_phe.py
                create_bin_phe_file(in_tsv = tsv, out_phe = phe, out_log = log, field_id = phe_values['field_id'],
                                    case     = phe_values['case'].replace(',',';').split(';'),
                                    control  = phe_values['control'].replace(',',';').split(';'),
                                    missing_is_control = all_ctrl if phe_values['control'] else True)
            else:
                create_qt_phe_file(in_tsv = tsv, out_phe = phe, out_log = log, field_id = phe_values['field_id'],
                                   order    = phe_values['qt_order'].replace(',',';').split(';'),
                                   exclude  = phe_values['exclude'].replace(',',';').split(';'))
            # annotate the phenotype
            make_phe_info([phe], 
                           info_dir, 
                          [phe_values['desc']],
                           phe_values['field_id'],
                           phe_values['table_id'],
                           phe_values['basket_id'],
                           '24983', 
                           os.path.basename(in_tsv))
            if os.path.exists(os.path.join(home_dir, "current", "phe/{}.phe".format(phe_name))):
                for folder, filetype in zip(["phe", "info", "logs"], ["phe", "info", "log"]):
                    os.remove(os.path.join(home_dir, "current", folder, phe_name + "." + filetype))
            for path, folder, filetype in zip([phe, info_file, log], ["phe", "info", "logs"], ["phe", "info", "log"]):
                os.symlink(path, os.path.join(home_dir, "current", folder, phe_name + "." + filetype))
    return

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=_README_
    )    
    parser.add_argument('--tsv', dest="input", required=True, nargs=1,
                            help='input table from phenotyping session')
    parser.add_argument('--no-header', dest="noheader", action='store_true',
                            help='flag if input tsv has no header')
    parser.add_argument('--name', dest="name", required=True, nargs=1, 
                            help='column in tsv corresponding to phe file name (GBE ID)')
    parser.add_argument('--desc', dest="desc", required=True, nargs=1, 
                            help='column in tsv corresponding to phenotype description (GBE NAME string)')
    parser.add_argument('--field', dest="field", required=True, nargs=1,
                            help='column in tsv corresponding to UK Biobank Field ID')
    parser.add_argument('--case', dest='case', required=True, nargs=1,
                            help='column in tsv corresponding to values for binary case definitions')
    parser.add_argument('--control', dest='control', required=True, nargs=1,
                            help='column in tsv corresponding to values for binary control definitions')
    parser.add_argument('--missing-is-control', dest="expand_control", action='store_true',
                            help='flag if missing values for binary traits should be defined as controls')
    parser.add_argument('--missing', dest='exclude', required=True, nargs=1,
                            help='column in tsv corresponding to QT values considered as missing data')
    parser.add_argument('--order', dest='order', required=True, nargs=1,
                            help='column in tsv corresponding to order of values (least to greatest) for QTs from categorical fields')
    parser.add_argument('--annotator', dest='annotator', required=True, nargs=1,
                            help='column in tsv corresponding to the annotator name')
    parser.add_argument('--annotator_date', dest='annotator_date', required=True, nargs=1,
                            help='column in tsv corresponding to the annotation date')
    parser.add_argument('--only-this-row', dest='onlyone', required=False, nargs=1,
                            help='(optional) flag to run only one (zero-indexed, not including the header) row the input tsv.') 
    args = parser.parse_args()
    # File gives the mapping between table, field, basket, and release date 
    ftbd = pd.read_csv('/oak/stanford/groups/mrivas/users/guhan/repos/ukbb-tools/05_gbe/field_table_basket_date.tsv', sep='\t', dtype=object)
    ftbd['Release_Date'] = pd.to_datetime(ftbd.Release_Date)
    # Retain only the latest release per field
    ftbd = ftbd.sort_values('Release_Date', ascending=False).groupby('Field_ID').first().reset_index()
    define_phenos(in_tsv = args.input[0],
           field_col = int(args.field[0]),
           name_col  = int(args.name[0]),
           case_col  = int(args.case[0]),
           ctrl_col  = int(args.control[0]),
           excl_col  = int(args.exclude[0]),
           qtfc_col  = int(args.order[0]),
           desc_col  = int(args.desc[0]),
           ann_col   = int(args.annotator[0]),
           ann_d_col = int(args.annotator_date[0]),
           ftbd      = ftbd,
           header    = not args.noheader,
           all_ctrl  = args.expand_control,
           make_this_one = args.onlyone)

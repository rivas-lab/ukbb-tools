#!/bin/python
import os,glob
import pandas as pd
from make_phe import *
from annotate_phe import make_phe_info

_README_="""
This script is designed to process table-defined phenotypes (from a Rivas lab computing session) into *.phe files usable for downstream analysis with PLINK and similar tools. 

To define all phenotypes from a table, identify the (zero-indexed!) columns specified below, and pass them to the script. If you would like to define only one phenotype from the table, you can do that by adding the flag --only-this-row, followed by the row (zero-indexed, such that the first non-header row is number 0) you'd like to run.

Author: Matthew Aguirre (SUNET: magu)
"""

def get_phe_definitions(in_tsv, header=True, special_row=None):
    # get phenotype definitions from an input table, with option to only extract one row
    phe_info = {}
    if make_this_one is not None:
        if isinstance(make_this_one, list):
            make_this_one = int(make_this_one[0])
        else:
            make_this_one = int(make_this_one)
    with open(in_tsv, 'r') as f:
        for i,line in enumerate(f):
            if header and i==0:
                if make_this_one is not None:
                    make_this_one += 1
                continue
            if make_this_one is not None and i != make_this_one:
                continue
            fields = line.rstrip().split('\t')
            if not fields[name_col]: continue
            # the if switch is in case empty cells get pruned from the TSV
            phe_info[ fields[name_col] ] = {'case':     fields[case_col] if case_col < len(fields) else '',
                                            'control':  fields[ctrl_col] if ctrl_col < len(fields) else '',
                                            'qt_order': fields[qtfc_col] if qtfc_col < len(fields) else '',
                                            'exclude':  fields[excl_col] if excl_col < len(fields) else '',
                                            'table_id': fields[table_col] if table_col < len(fields) else '',
                                            'field_id': fields[field_col] if field_col < len(fields) else '',
                                            'desc':     fields[desc_col] if desc_col < len(fields) else ''}
    return phe_info


def define_phenos(in_tsv, table_col, field_col, name_col, case_col, ctrl_col, 
                  excl_col, qtfc_col, desc_col, new_table=True, table_id=None, 
                  header=True, all_ctrl=False, make_this_one=None):
    home_out_dir='/oak/stanford/groups/mrivas/private_data/ukbb/24983/phenotypedata/'
    home_in_dir ='/oak/stanford/groups/mrivas/private_data/ukbb/24983/phenotypedata/download/'
    # iterate over phenotypes in the input file
    for phe_name, phe_values in get_phe_definitions(in_tsv, header, make_this_one).items():
        print(phe_info)
        print(phe_name, phe_values)
        # select tab file to use from handler arguments
        if new_table:
            tsv = map(glob.glob, [os.path.join(root,'newest/ukb*.tab') for root,dirs,files in os.walk(home_in_dir) if phe_values['table_id'] in dirs])[0][0]
            table_id = os.path.splitext(os.path.basename(tsv))[0].replace('ukb','')
        else:
            table_id = phe_values['table_id'] if table_id is None else table_id
            tab_f = 'ukb{}.tab'.format(table_id)
            # this will throw an indexing error if a bad table is supplied
            tsv = [os.path.join(root,tab_f) for root,dirs,files in os.walk(home_in_dir) if tab_f in files][0]
        # get phenotype name
        basket_id = os.path.basename(os.path.dirname(os.path.dirname(tsv)))
        phe = os.path.join(home_out_dir, basket_id, table_id, phe_name+'.phe')
        print(tsv,phe)
        log = os.path.join(os.path.dirname(phe), "logs/{0}.log".format(phe_name))
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
                       os.path.join(os.path.dirname(phe), "info"), 
                      [phe_values['desc']],
                       phe_values['field_id'],
                       phe_values['table_id'],
                       basket_id,
                       '24983', 
                       os.path.basename(in_tsv))
                       
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
    parser.add_argument('--table', dest="table", required=True, nargs=1,
                            help='column in tsv corresponding to UK Biobank Table ID')
    parser.add_argument('--table-id', dest="table_id", required=False, nargs=1, default=[None],
                            help='ID of UK Biobank Table ID to use (overrides table option)')
    parser.add_argument('--table-newest', dest="new_table", action='store_true',
                            help='Flag to use newest table corresponding to UK Biobank basket (overrides table and table-id options)')
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
    parser.add_argument('--only-this-row', dest='onlyone', required=False, nargs=1,
                            help='(optional) flag to run only one (zero-indexed, not including the header) row the input tsv.') 
    args = parser.parse_args()
    print(args)
    # lol i hope this works
    define_phenos(in_tsv    = args.input[0],
           table_col = int(args.table[0]),
           field_col = int(args.field[0]),
           name_col  = int(args.name[0]),
           case_col  = int(args.case[0]),
           ctrl_col  = int(args.control[0]),
           excl_col  = int(args.exclude[0]),
           qtfc_col  = int(args.order[0]),
           desc_col  = int(args.desc[0]),
           new_table = args.new_table,
           table_id  = args.table_id[0],
           header    = not args.noheader,
           all_ctrl  = args.expand_control,
           make_this_one = args.onlyone)

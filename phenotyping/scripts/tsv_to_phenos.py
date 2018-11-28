#!/bin/python
import os
from make_phe import *

_README_="""
This script is designed to process table-defined phenotypes (from a Rivas lab computing session) into *.phe files usable for downstream analysis with PLINK and similar tools. 

To define all phenotypes from a table, identify the (zero-indexed!) columns specified below, and pass them to the script. If you would like to define only one phenotype from the table, you can do that by adding the flag --only-this-row, followed by the row (zero-indexed, such that the first non-header row is number 0) you'd like to run.

Author: Matthew Aguirre (SUNET: magu)
"""

def table_to_basket(table):
    basket_to_table = {'9796':['9797', '21732', '24611'],
                       '10136':['10137', '21731'],
                       '10483':['10484', '21733', '24613'],
                       '11139':['11140', '21734', '24614'],
                       '2000269':['21730', '24615'],
                       '2001702':['25279']
                      }
    for basket, tables in basket_to_table.items():
        if table in tables:
            return basket
    raise ValueError("Table {} not found in basket_to_table!".format(table))

def make_table(in_tsv, table_col, field_col, name_col, case_col, ctrl_col, 
               excl_col, qtfc_col, header=True, all_ctrl=False, make_this_one=None):
    home_out_dir='/oak/stanford/groups/mrivas/dev-ukbb-tools/phenotypes/'
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
                                            'field_id': fields[field_col] if field_col < len(fields) else ''}
    # this info should be logged somewhere
    for phe_name, phe_values in phe_info.items():
        print(phe_name, phe_values)
        if phe_values['case']: # assume binary if we have a case definition
            # this and create_qt_phe_file below are implemented in make_phe.py
            create_bin_phe_file(in_tsv   = os.path.join("/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/download",
                                                       table_to_basket(phe_values['table_id']),
                                                       phe_values['table_id'], 
                                                       "ukb{}.tab".format(phe_values['table_id'])),
                                out_phe  = os.path.join(home_out_dir, phe_values['table_id'], "{0}.phe".format(phe_name)),
                                out_log  = os.path.join(home_out_dir, phe_values['table_id'], "logs", "{0}.log".format(phe_name)),
                                field_id = phe_values['field_id'],
                                case     = phe_values['case'].replace(',',';').split(';'),
                                control  = phe_values['control'].replace(',',';').split(';'),
                                missing_is_control = all_ctrl)
        else: # assume qt
            create_qt_phe_file(in_tsv   = os.path.join("/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/download",
                                                       table_to_basket(phe_values['table_id']), 
                                                       phe_values['table_id'], 
                                                       "ukb{}.tab".format(phe_values['table_id'])),
                               out_phe  = os.path.join(home_out_dir, phe_values['table_id'], "{0}.phe".format(phe_name)),
                               out_log  = os.path.join(home_out_dir, phe_values['table_id'], "logs", "{0}.log".format(phe_name)),
                               field_id = phe_values['field_id'],
                               order    = phe_values['qt_order'].replace(',',';').split(';'),
                               exclude  = phe_values['exclude'].replace(',',';').split(';'))

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
    parser.add_argument('--field', dest="field", required=True, nargs=1,
                            help='column in tsv corresponding to UK Biobank Field ID')
    parser.add_argument('--table', dest="table", required=True, nargs=1,
                            help='column in tsv corresponding to UK Biobank Table ID')
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
    make_table(in_tsv    = args.input[0],
           table_col = int(args.table[0]),
           field_col = int(args.field[0]),
           name_col  = int(args.name[0]),
           case_col  = int(args.case[0]),
           ctrl_col  = int(args.control[0]),
           excl_col  = int(args.exclude[0]),
           qtfc_col  = int(args.order[0]),
           header    = not args.noheader,
           all_ctrl  = args.expand_control,
           make_this_one = args.onlyone)

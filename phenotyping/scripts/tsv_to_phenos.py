#!/bin/python
import os,glob
import pandas as pd
from make_phe import *

_README_="""
This script is designed to process table-defined phenotypes (from a Rivas lab computing session) into *.phe files usable for downstream analysis with PLINK and similar tools. 

To define all phenotypes from a table, identify the (zero-indexed!) columns specified below, and pass them to the script. If you would like to define only one phenotype from the table, you can do that by adding the flag --only-this-row, followed by the row (zero-indexed, such that the first non-header row is number 0) you'd like to run.

Author: Matthew Aguirre (SUNET: magu)
"""

def make_table(in_tsv, table_col, field_col, name_col, case_col, ctrl_col, 
               excl_col, qtfc_col, desc_col, 
               header=True, all_ctrl=False, make_this_one=None):
    home_out_dir='/oak/stanford/groups/mrivas/dev-ukbb-tools/phenotypes/'
    home_in_dir ='/oak/stanford/groups/mrivas/private_data/ukbb/24983/phenotypedata/download/'
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
    # this info should be logged somewhere
    for phe_name, phe_values in phe_info.items():
        print(phe_name, phe_values)
        # these are the same regardless of the nature of the phenotype
        tab = phe_values['table_id'] # for brevity below
        tsv = map(glob.glob, [os.path.join(root,'newest/ukb*.tab') for root,dirs,files in os.walk(home_in_dir) if tab in dirs])[0][0]
        # home_out_dir/basketID/gbeID.phe
        phe = os.path.join(home_out_dir, os.path.basename(os.path.dirname(os.path.dirname(tsv))), '{0}.phe'.format(phe_name))
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
        update_phe_info(os.path.splitext(log)[0]+'.info', phe_values['desc'], in_tsv) 
    return


def update_phe_info(info_f, phe_desc, source_table):
    info = pd.read_table(info_f, index_col=0)
    info.insert(0,  "GBE_NAME", [phe_desc])
    info.insert(10, "SOURCE", [os.path.basename(source_table)])
    info.to_csv(info_f, sep="\t")
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
           desc_col  = int(args.desc[0]),
           header    = not args.noheader,
           all_ctrl  = args.expand_control,
           make_this_one = args.onlyone)

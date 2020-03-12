#!/bin/python
#coding=utf-8
import os
import numpy as np
import pandas as pd
import glob
from make_phe import *
from annotate_phe import make_phe_info
from showcase_and_list_to_tsv import join_and_add_cols

_README_ = """
A script to tell which fields in a UK Biobank table (.tab) file are new or updated.

Actions:
1. If table contains new data fields, make table for computing session
2. If table contains updated data fields:
    a. Make updated .phe files (and updates link to most recent version — this is TODO)
    b. Run gwas on said .phe files
    c. Optionally, suppress either of the above with flags described below.

Authors: Matthew Aguirre and Guhan Venkataraman(SUNETs: magu and guhan)
Updated: 2020/01/29
"""
new_col_order = ['Annotator', 'Annotation date', 'Name', 'GBE ID', 'Field', 'FieldID',
        'QT_total_num', 'BIN_total_num', 'QT_index', 'BIN_index', 'coding_exclude', 'coding_QT',
        'coding_binary_case', 'coding_binary_control', 'Participants', 'Stability', 'ValueType',
        'Units', 'Strata', 'Sexed', 'Instances', 'Array', 'Coding', 'Link']


def create_tsv(new_f):
    complete_new_header_df = pd.read_table(new_f, sep="\t", index_col='f.eid', nrows=1)
    new_fields = set([s.split('.')[1] for s in complete_new_header_df.columns])
    # make table for computing session with showcase_and_list_to_tsv 
    if len(new_fields) != 0:
        out_file = '../tables/ukb_' + str(pd.datetime.today().date()).replace('-','') + '.tsv'
        out_df = join_and_add_cols(map(int, [x for x in iter(new_fields)]))
        # Add in all previous annotations
        tsvs = glob.glob('/oak/stanford/groups/mrivas/users/guhan/repos/ukbb-tools/02_phenotyping/tables/*.tsv')
        tsvs = [tsv for tsv in tsvs if ((not 'params' in tsv) and (not 'priority' in tsv))]
        for tsv in tsvs:
            out_df.update(pd.read_table(tsv, dtype=object))
        out_df.to_csv(out_file, sep='\t', index=False)
        print("New .tsv made: " + out_file)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=_README_
    )
    parser.add_argument('--basket', dest='basket', required=True,
                        help='input basket id')
    parser.add_argument('--table', dest='table', required=True,
                        help='input table id/path')
    args = parser.parse_args()
    in_b = args.basket
    if os.path.exists(args.table):
        new_f = os.path.abspath(args.table)
    else:
        new_f = os.path.abspath(os.path.join("/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/",
                             in_b, args.table, "download", 'ukb'+args.table+'.tab'))
        if not os.path.exists(new_f):
            raise OSError("Could not find input file {0}".format(new_f))
    # 2. determine which table to compare input to 
    print("Creating .tsv for new table for basket {0}:\n {1}...\n".format(in_b, new_f))
    create_tsv(new_f)

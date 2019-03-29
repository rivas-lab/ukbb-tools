#!/bin/python
import os
import pandas as pd
from make_phe import *

# actions:
#   1. if new data, output new fields and make table for computing session
#   2. if not new data:
#       a. figure out which fields (if any) are new, apply #1
#       b. among old fields, figure out which ones have been updated (higher counts)
#           i. optionally, update these in the gbe pipeline (analysis + gwas)


_README_ = """
A script to tell which fields in a UK Biobank table (.tab) file are new or updated.

Outputs:
1. A list of new fields, and spreadsheet* for phenotype definition computing session.
2. New phenotype files* and summary stats* from array gwas using updated data.

Starred items can be suppressed. See flags for more info.

Author: Matthew Aguirre (SUNET: magu)
"""

def find_new_data(new_f, old_f, make_table):
    # get new fields
    new_df = pd.read_table(new_f, sep="\t", index_col='f.eid', nrows=5)
    old_df = pd.read_table(old_f, sep="\t", index_col='f.eid', nrows=5)
    new_df_col = set([s.split('.')[1] for s in new_df.columns])
    old_df_col = set([s.split('.')[1] for s in old_df.columns])
    new_fields = {i for i in iter(new_df_col) if i not in old_df_col}
    old_fields = {i for i in iter(new_df_col) if i in old_df_col}
    print("New fields:\n" + "\n".join(iter(new_fields)))
    print("Updated fields:")
    # make table for computing session with 
    if len(new_fields) != 0 and make_table:
        from showcase_and_list_to_tsv import join_and_add_cols
        out_file = '../tables/ukb_' + str(pd.datetime.today().date()).replace('-','') + '.tsv'
        join_and_add_cols(new_fields).to_csv(out_file, sep='\t', index=False)
    # determine which fields got updated
    updated_fields = []
    for field in iter(old_fields):
        # load up columns for this field
        new_df = pd.read_table(new_f, sep="\t", index_col='f.eid',
                               usecols=lambda s: s=='f.eid' or s.split('.')[1]==field
                               ).sort_index()
        old_df = pd.read_table(old_f, sep="\t", index_col='f.eid', 
                               usecols=lambda s: s=='f.eid' or s.split('.')[1]==field
                               ).sort_index()
        # subset to same set of individuals, in case some are redacted in new data
        shared_inds = new_df.index.intersection(old_df.index)
        # compare values
        if not new_df.loc[shared_inds,:].equals(old_df.loc[shared_inds,:]):
            print(field)
            updated_fields.append(field)
    return updated_fields

def update_phenos(fields, ukb_tab, table_id, basket_id):
    # iterate over all references
    import glob
    phe_data_root = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/'
    table_info = pd.read_table('../tables/gbe_sh_input_params.tsv', index_col=0)
    for gbe_table in table_info.index:
        name_col = table_info.loc[gbe_table,'nameCol (GBE ID)']
        field_col= table_info.loc[gbe_table,'fieldCol (field_ID)']
        case_col = table_info.loc[gbe_table,'caseCol (coding_binary_case)']
        ctrl_col = table_info.loc[gbe_table,'ctrlCol (coding_binary_control)']
        excl_col = table_info.loc[gbe_table,'exclCol (coding_exclude)']
        order_col= table_info.loc[gbe_table,'orderCol (coding_QT)']
        desc_col = table_info.loc[gbe_table, 'descCol']
        phe_defs = pd.read_table('../tables/'+gbe_table)
        for phe in phe_defs.loc[phe_defs.iloc[field_col] in fields,:].rows:
            phe_file = phe_defs['nameCol (GBE ID)'] + '.phe'
            phe_log = phe_defs['nameCol (GBE ID)'] + '.log'
            if not phe[case_col].isnull():
                create_bin_phe_file(in_tsv = ukb_tab,
                                    out_phe= os.path.join(phe_data_root,basket_id,table_id,phe_file),
                                    out_log= os.path.join(phe_data_root,basket_id,table_id,'logs',phe_logs),
                                    field_id = phe[field_col],
                                    case   = phe[case_col],
                                    control= phe[ctrl_col],
                                    missing_is_control=False
                                    )
            else:
                create_qt_phe_file(in_tsv = ukb_tab,
                                   out_phe= os.path.join(phe_data_root,basket_id,table_id,phe_file),
                                   out_log= os.path.join(phe_data_root,basket_id,table_id,'logs',phe_logs),
                                   field_id = phe[field_col],
                                   exclude= phe[excl_col],
                                   order=phe[order_col]
                                   )
    return

def update_summary_stats(fields):
    return 'todo'

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
    parser.add_argument('--old-table', dest='old_table', required=False, default=None,
                        help='old table id/path')   
    parser.add_argument('--no-update', dest='no_update', action='store_true', default=False,
                        help='flag to not automagically update phenotypes or run array gwas')
    parser.add_argument('--no-gwas', dest='no_gwas', action='store_true', default=False,
                        help='flag to not automagically run array gwas on updated phenotypes')
    parser.add_argument('--no-spreadsheet', dest='no_table', action='store_true', default=False,
                        help='flag to not automagically create a table from new phenotypes for computing session')
    args = parser.parse_args()
    print(args)
     
    # 1. determine if table is new 
    in_b = args.basket
    if os.path.exists(args.table):
        new_f = args.table
    else:
        new_f = os.path.join("/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/download",
                             in_b, args.table, 'ukb'+args.table+'.tab')
        if not os.path.exists(new_f):
            raise OSError("Could not find input file {0}".format(new_f))
    # 2. determine which table to compare input to 
    if os.path.exists(args.old_table):
        old_f = args.old_table
    elif args.old_table is not None:
        old_f = os.path.join("/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/download", in_b, args.old_table, "raw.tsv")
    else:
        old_f = os.path.join("/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/download", in_b, "newest", "raw.tsv")
    print("Analyzing new table for basket {0}:\n {1}...\n".format(in_b, new_f))
    
    if not os.path.exists(old_f):
        print("could not find existing table for basket {0}, assuming this is a new basket!\n".format(in_b))
      
    # 3. diff files to get new (printed/written to file) and updated fields (stored as variable)
    fields = find_new_data(new_f=new_f, old_f=old_f, make_table=not args.no_table)
    print(fields)
    # 4. update phenotypes
    if not args.no_update:
        update_phenos(fields=fields, ukb_tab=new_f, table_id=, basket_id=args.basket)

     

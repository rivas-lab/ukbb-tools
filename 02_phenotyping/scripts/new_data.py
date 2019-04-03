#!/bin/python
import os
import pandas as pd
from make_phe import *
from annotate_phe import make_phe_info
import numpy as np

_README_ = """
A script to tell which fields in a UK Biobank table (.tab) file are new or updated.

Actions:
1. If table contains new data fields, make table for computing session
2. If table contains updated data fields:
    a. Make updated .phe files (and updates link to most recent version — this is TODO)
    b. Run gwas on said .phe files
    c. Optionally, suppress either of the above with flags described below.

Author: Matthew Aguirre (SUNET: magu)
"""

# TODO: add flag to update master phenotype file

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
            return [field]
            updated_fields.append(field)
    return updated_fields


def update_phenos(fields, ukb_tab, table_id, basket_id):
    # get data and requisite info
    phe_data_root = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/'
    table_info = pd.read_table('../tables/gbe_sh_input_params.tsv', index_col=0, dtype=str)
    paths_to_phenos = []
    phe_data_root = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/'
    # call tsv_to_phenos with --only-this-row to get the phenotypedata
    for gbe_table in table_info.index:
        phe_defs = pd.read_table('../tables/'+gbe_table, dtype=str).fillna('')
        name_col = int(table_info.loc[gbe_table,'nameCol (GBE ID)'])
        field_col= int(table_info.loc[gbe_table,'fieldCol (field_ID)'])
        for ix,phe in phe_defs.loc[phe_defs.iloc[:,field_col].isin(fields),:].iterrows():
            os.system(' '.join(('python tsv_to_phenos.py', 
                           '--tsv', '../tables/'+gbe_table,
                           '--table', table_info.loc[gbe_table, 'tableCol (table_ID)'],
                           '--table-id', table_id,
                           '--name', table_info.loc[gbe_table,'nameCol (GBE ID)'],
                           '--desc', table_info.loc[gbe_table, 'descCol'],
                           '--field', table_info.loc[gbe_table,'fieldCol (field_ID)'],
                           '--case', table_info.loc[gbe_table,'caseCol (coding_binary_case)'],
                           '--control', table_info.loc[gbe_table,'ctrlCol (coding_binary_control)'],
                           '--order', table_info.loc[gbe_table,'orderCol (coding_QT)'],
                           '--missing', table_info.loc[gbe_table,'exclCol (coding_exclude)'],
                           '--only-this-row', str(ix))))
            phe_file = os.path.join(phe_data_root,basket_id,table_id,phe[name_col] + '.phe')
            paths_to_phenos.append(phe_file)
    return paths_to_phenos

 
def update_summary_stats(phe_files):
    for f in phe_files:
        os.system(" ".join(["python ../../04_gwas/gwas.py --run-array",
                                   "--pheno", f, "--population white_british",
                                   "--log-dir", os.path.join(os.path.dirname(f).replace('phenotypedata','cal/gwas'), 'logs'),
                                   "--out", os.path.dirname(f).replace('phenotypedata','cal/gwas')]))
    return

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
        phe_files = update_phenos(fields=fields, 
                                  ukb_tab=new_f, 
                                  table_id=os.path.splitext(os.path.basename(new_f))[0].replace('ukb',''),
                                  basket_id=args.basket)
    else:
        phe_files = []
    # 5. update gwas
    if not args.no_gwas and phe_files:
        update_summary_stats(phe_files)
    

#!/bin/python
# coding=utf-8
import os
import numpy as np
import pandas as pd
import glob
from make_phe import *
from annotate_phe import make_phe_info
from showcase_and_list_to_tsv import join_and_add_cols

_README_ = """

A script that generates a new master annotation tsv.

Actions:
1. Gets a tsv mapping fields to tables, baskets, and release dates.
2. Gets most recent version of each field.
3. Gets annotations or empty fields generated thus far.

Authors: Matthew Aguirre and Guhan Venkataraman (SUNETs: magu and guhan)
Updated: 2020/03/16

"""

def create_tsv(new_col_order):
    ftbd = pd.read_csv('../tables/field_table_basket_date.tsv', sep='\t')
    all_fields = set(ftbd['FieldID'])
    # make table for computing session with showcase_and_list_to_tsv
    out_file = (
        "../tables/annotations/ukb_" + str(pd.datetime.today().date()).replace("-", "") + ".tsv"
    )
    out_df = join_and_add_cols(map(int, [x for x in iter(all_fields)])).astype(
        object
    )
    # Gather all previous annotations
    tsvs = glob.glob(
        "/oak/stanford/groups/mrivas/users/guhan/repos/ukbb-tools/02_phenotyping/tables/annotations/*.tsv"
    )
    tsvs = [
        pd.read_csv(tsv, sep="\t", dtype={"FieldID": int})
        for tsv in tsvs
    ]
    prev_annots = pd.concat(tsvs, sort=True)[new_col_order]
    prev_annots = prev_annots.merge(out_df, how="left", on=["FieldID"])
    # Get annotations from the new tsv if not in previous table for some reason
    for col in new_col_order:
        if col != "FieldID":
            prev_annots[col] = prev_annots.apply(
                lambda row: row[col + "_y"]
                if pd.isna(row[col + "_x"])
                else row[col + "_x"],
                axis=1,
            )
    prev_annots = prev_annots[new_col_order]
    print(
        "Number of distinct previous annotations prior to dropping duplicates: "
        + str(len(set(prev_annots[prev_annots["GBE ID"].notna()]["GBE ID"])))
    )
    print(
        "Number of NA fields from before: "
        + str(len(set(prev_annots[prev_annots["GBE ID"].isna()]["FieldID"])))
    )
    # Get rid of duplicate previous annotations based on GBE ID/FieldID
    prev_annots = prev_annots.sort_values(new_col_order).drop_duplicates(
        ["GBE ID", "FieldID"]
    )
    print(
        "Number of distinct previous annotations after dropping duplicates: "
        + str(len(set(prev_annots[prev_annots["GBE ID"].notna()]["GBE ID"])))
    )
    print(
        "Number of NA fields from before after dropping duplicates: "
        + str(len(set(prev_annots[prev_annots["GBE ID"].isna()]["FieldID"])))
    )
    annotated = set(prev_annots["FieldID"])
    # Annotate them back in
    out_df = out_df[~out_df["FieldID"].isin(annotated)]
    out_df = pd.concat([out_df, prev_annots])
    out_df[["FieldID", "QT_total_num", "BIN_total_num", "QT_index", "BIN_index", "Participants", "Instances", "Array", "Coding"]] = out_df[["FieldID", "QT_total_num", "BIN_total_num", "QT_index", "BIN_index", "Participants", "Instances", "Array", "Coding"]].astype(float)
    out_df[["FieldID", "QT_total_num", "BIN_total_num", "QT_index", "BIN_index", "Participants", "Instances", "Array", "Coding"]] = out_df[["FieldID", "QT_total_num", "BIN_total_num", "QT_index", "BIN_index", "Participants", "Instances", "Array", "Coding"]].astype("Int64")
    out_df = out_df.drop_duplicates()
    out_df.to_csv(out_file, sep="\t", index=False)
    print("New .tsv made: " + out_file + ".")
    if os.path.exists("../tables/ukb_annotations.tsv"):
        os.remove("../tables/ukb_annotations.tsv")
    os.symlink(out_file, "../tables/ukb_annotations.tsv")
    print("Symlink for ../tables/ukb_annotations.tsv updated to " + out_file + ".")


if __name__ == "__main__":
    import argparse
    print("Updating ../tables/ukb_annotations.tsv...")
    new_col_order = [
        "Annotator",
        "Annotation date",
        "Name",
        "GBE ID",
        "Field",
        "FieldID",
        "QT_total_num",
        "BIN_total_num",
        "QT_index",
        "BIN_index",
        "coding_exclude",
        "coding_QT",
        "coding_binary_case",
        "coding_binary_control",
        "Participants",
        "Stability",
        "ValueType",
        "Units",
        "Strata",
        "Sexed",
        "Instances",
        "Array",
        "Coding",
        "Link",
    ]
    create_tsv(new_col_order)

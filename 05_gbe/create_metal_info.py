#!/bin/python

import pandas as pd
import os

pops = ["white_british", "non_british_white", "african", "e_asian", "s_asian", "related", "others"]
dfs = []

for pop in pops:
    df = pd.read_table('array-combined/icdinfo.array.' + pop + '.tsv')
    dfs.append(df)

ns = []
for gbe_id in list(df['GBE_ID']):
    # find out which populations to add
    log = "/oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current/metal/" + gbe_id + ".metal.info.txt"
    if os.path.exists(log):
        with open(log) as f:
            lines = f.read().splitlines()
        input_files = [line for line in lines if (("oak" in line) and ("private_data" not in line))]
        # add them up
        populations = [input_file.split("/")[-2] for input_file in input_files]
        print(gbe_id)
        print(populations)
        n = 0
        for df, pop in zip(dfs, pops):
            if pop in populations:
                n += df[df['GBE_ID'] == gbe_id]['N'].values[0]
        ns.append(str(n))
    else:
        ns.append("")

df = df.drop(columns="N")
df["N"] = ns
df["N"] = df["N"]
df = df[["GBE_category", "GBE_ID", "N", "GBE_NAME", "GBE_short_name", "GBE_short_name_len"]]
df.to_csv('array-combined/icdinfo.array.metal.tsv', sep='\t', index=False)

for pop in pops:
    df = pd.read_table('exome/200k/icdinfo.exome.' + pop + '.tsv')
    dfs.append(df)

ns = []
for gbe_id in list(df['GBE_ID']):
    # find out which populations to add
    log = "/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/current/metal/ukb24983_exomeOQFE." + gbe_id + ".metal.info.txt"
    if os.path.exists(log):
        with open(log) as f:
            lines = f.read().splitlines()
        input_files = [line for line in lines if (("scratch" in line) and ("metal" not in line))]
        # add them up
        populations = [input_file.split("/")[-2] for input_file in input_files]
        print(gbe_id)
        print(populations)
        n = 0
        for df, pop in zip(dfs, pops):
            if pop in populations:
                n += df[df['GBE_ID'] == gbe_id]['N'].values[0]
        ns.append(str(n))
    else:
        ns.append("")

df = df.drop(columns="N")
df["N"] = ns
df["N"] = df["N"]
df = df[["GBE_category", "GBE_ID", "N", "GBE_NAME", "GBE_short_name", "GBE_short_name_len"]]
df.to_csv('exome/200k/icdinfo.exome.metal.tsv', sep='\t', index=False)

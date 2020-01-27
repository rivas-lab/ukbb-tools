import numpy as np
import pandas as pd
import pickle
import os

sumstat_paths = pd.read_table('sumstat_paths.tsv')
phenos = list(sumstat_paths['PHEN'])

ptv = {analysis: [] for analysis in ["wb", "ma"]}
pav = {analysis: [] for analysis in ["wb", "ma"]}

included = 0

for i, pheno in enumerate(phenos):
    if i % 100 == 0:
        print("Done " + str(i) + " phenotypes out of 2537")
    stream = os.popen('find /oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_rv_array -name "*_' + pheno + '_*"')
    wb_path = stream.read().strip()
    stream = os.popen('find /oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_rv_ma_array -name "*_' + pheno + '_*"')
    ma_path = stream.read().strip()
    if (wb_path != "") and (ma_path != ""):
        wb, ma = pd.read_table(wb_path), pd.read_table(ma_path)
        if len(wb) > 1 and len(ma) > 1:
            wb_columns_to_rename, ma_columns_to_rename = wb.columns[1:], ma.columns[1:]
            renamed_columns_wb, renamed_columns_ma = [(col + "_wb") for col in wb_columns_to_rename], [(col + "_ma") for col in ma_columns_to_rename]
            wb.rename(columns=dict(zip(wb_columns_to_rename, renamed_columns_wb)), inplace=True)
            ma.rename(columns=dict(zip(ma_columns_to_rename, renamed_columns_ma)), inplace=True)
            merged = wb.merge(ma, on="gene")
            for var_type, bflist in zip(["pav", "ptv"], [ptv, pav]):
                wb_col = "log_10_BF_study_similar_var_independent_sigma_m_var_" + var_type + "_wb"
                ma_col = "log_10_BF_study_similar_var_independent_sigma_m_var_" + var_type + "_ma"
                if (ma_col in merged.columns) and (wb_col in merged.columns):
                    included += 1
                    tmp = merged[(merged[wb_col] >= 5) | (merged[ma_col] >= 5)]
                    bflist["wb"].extend(list(tmp[wb_col]))
                    bflist["ma"].extend(list(tmp[ma_col]))

for var_type, bflist in zip(["pav", "ptv"], [ptv, pav]):
    f = open(var_type + ".pkl","wb")
    pickle.dump(bflist, f)
    f.close()

print(str(included/2) + " phenotypes included in dictionaries.")


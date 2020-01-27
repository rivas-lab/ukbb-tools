import numpy as np
import pandas as pd

phen_info = pd.read_table("../05_gbe/phenotype_info.tsv")
phenos = list(phen_info[(phen_info["N"] >= 10) & ~(phen_info["#GBE_ID"].str.startswith("MED"))]["#GBE_ID"])

data = []

print("loading in data")
for i in range(len(phenos)):
    pearson = np.load("R_phen_rows/" + str(i + 1) + "_pearson.npy")
    pearson_p = np.load("R_phen_rows/" + str(i + 1) + "_pearson_p.npy")
    pearson_nohla = np.load("R_phen_rows/" + str(i + 1) + "_pearson_nohla.npy")
    pearson_p_nohla = np.load("R_phen_rows/" + str(i + 1) + "_pearson_p_nohla.npy")
    spearman = np.load("R_phen_rows/" + str(i + 1) + "_spearman.npy")
    spearman_p = np.load("R_phen_rows/" + str(i + 1) + "_spearman_p.npy")
    spearman_nohla = np.load("R_phen_rows/" + str(i + 1) + "_spearman_nohla.npy")
    spearman_p_nohla = np.load("R_phen_rows/" + str(i + 1) + "_spearman_p_nohla.npy")
    for j in range(len(phenos)):
        if i < j:
            data.append([phenos[i], phenos[j], pearson[j], pearson_p[j], pearson_nohla[j], pearson_p_nohla[j], spearman[j], spearman_p[j], spearman_nohla[j], spearman_p_nohla[j]])
        elif i == j:
            data.append([phenos[i], phenos[j], 1, 0, 1, 0, 1, 0, 1, 0])
                
            
df = pd.DataFrame(data, columns=['PHEN_1', 'PHEN_2', 'PEARSON', 'PEARSON_P', 'PEARSON_NOHLA', 'PEARSON_P_NOHLA', 'SPEARMAN', 'SPEARMAN_P', 'SPEARMAN_NOHLA', 'SPEARMAN_P_NOHLA'])

#df = df.dropna()

short_names = pd.read_table('../05_gbe/icdinfo.white_british.shortname.txt', header=None, names=["PHEN", "N1", "longname", "N2", "N3", "BIN", "shortname", "case"])

print("merging on short names")
df = df.merge(short_names, left_on="PHEN_1", right_on="PHEN", how="left")
df = df[['PHEN_1', 'shortname', 'PHEN_2', 'PEARSON', 'PEARSON_P', 'PEARSON_NOHLA', 'PEARSON_P_NOHLA', 'SPEARMAN', 'SPEARMAN_P', 'SPEARMAN_NOHLA', 'SPEARMAN_P_NOHLA']]
df.columns = ['PHEN_1', 'shortname_1', 'PHEN_2', 'PEARSON', 'PEARSON_P', 'PEARSON_NOHLA', 'PEARSON_P_NOHLA', 'SPEARMAN', 'SPEARMAN_P', 'SPEARMAN_NOHLA', 'SPEARMAN_P_NOHLA']

df = df.merge(short_names, left_on="PHEN_2", right_on="PHEN", how="left")
df = df[['PHEN_1', 'shortname_1', 'PHEN_2', 'shortname', 'PEARSON', 'PEARSON_P', 'PEARSON_NOHLA', 'PEARSON_P_NOHLA', 'SPEARMAN', 'SPEARMAN_P', 'SPEARMAN_NOHLA', 'SPEARMAN_P_NOHLA']]
df.columns = ['PHEN_1', 'shortname_1', 'PHEN_2', 'shortname_2', 'PEARSON', 'PEARSON_P', 'PEARSON_NOHLA', 'PEARSON_P_NOHLA', 'SPEARMAN', 'SPEARMAN_P', 'SPEARMAN_NOHLA', 'SPEARMAN_P_NOHLA']

df.to_csv('corrs_names.tsv', sep='\t', index=False)

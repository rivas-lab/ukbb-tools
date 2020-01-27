import numpy as np
import pandas as pd

phen_info = pd.read_table("../05_gbe/phenotype_info.tsv")
phenos = list(phen_info[(phen_info["N"] >= 10) & ~(phen_info["#GBE_ID"].str.startswith("MED"))]["#GBE_ID"])

data = []

print("loading in data")
for i in range(len(phenos)):
    betalen = np.load("R_phen_rows/" + str(i + 1) + "_betalen.npy")
    betalen_nohla = np.load("R_phen_rows/" + str(i + 1) + "_betalen_nohla.npy")
    for j in range(len(phenos)):
        if i <= j:
            data.append([phenos[i], phenos[j], betalen[j], betalen_nohla[j]])
            
df = pd.DataFrame(data, columns=['PHEN_1', 'PHEN_2', 'BETALEN', 'BETALEN_NOHLA'])

short_names = pd.read_table('../05_gbe/icdinfo.white_british.shortname.txt', header=None, names=["PHEN", "N1", "longname", "N2", "N3", "BIN", "shortname", "case"])

print("merging on short names")
df = df.merge(short_names, left_on="PHEN_1", right_on="PHEN", how="left")
df = df[['PHEN_1', 'shortname', 'PHEN_2', 'BETALEN', 'BETALEN_NOHLA']]
df.columns = ['PHEN_1', 'shortname_1', 'PHEN_2', 'BETALEN', 'BETALEN_NOHLA']

df = df.merge(short_names, left_on="PHEN_2", right_on="PHEN", how="left")
df = df[['PHEN_1', 'shortname_1', 'PHEN_2', 'shortname', 'BETALEN', 'BETALEN_NOHLA']]
df.columns = ['PHEN_1', 'shortname_1', 'PHEN_2', 'shortname_2', 'BETALEN', 'BETALEN_NOHLA']

df.to_csv('betalens_names.tsv', sep='\t', index=False)

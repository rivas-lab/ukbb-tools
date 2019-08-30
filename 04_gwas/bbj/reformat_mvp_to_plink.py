from __future__ import print_function
from __future__ import division 
import scipy.stats as st
import os, sys
import gzip 
import numpy as np
fin = open(sys.argv[1],'r').readlines()
bq = sys.argv[2]

fout = gzip.open(os.path.join(sys.arv[1], '.plink.gz'), 'wb')

for line in fin[0:]:
    line = line.rstrip()
    line = line.split('\t')
    if len(line) < 22:
        continue
    if line[0] == "ID":
        if bq == "BIN":
            #CHROM  POS     ID      REF     ALT     A1      FIRTH?  TEST    OBS_CT  OR      SE      Z_STAT  P
            print('#CHROM\tPOS\tID\tREF\tALT\tA1\tFIRTH?\tTEST\tOBS_CT\tOR\tSE\tZ_STAT\tP', file = fout)
        elif bq == "QT":
            #CHROM  POS     ID      REF     ALT     A1      TEST    OBS_CT  BETA    SE      T_STAT  P
            print('#CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP', file = fout)
    else:
        if len(line[17]) < 1:
            line[17] = float(line[16])/st.norm.ppf(1-float(line[3]))
            if bq == "BIN":
                print(line[6],line[7],line[2], line[11], line[12], line[19],"N",  "ADD", line[20], np.exp(line[16]), line[17], float(line[16])/float(line[17]), line[3], sep = '\t',  file = fout)
            elif bq == "QT":
                print(line[6],line[7],line[2], line[11], line[12], line[19], "ADD", line[20], line[16], line[17], float(line[16])/float(line[17]), line[3], sep = '\t',  file = fout)
			
	
fout.close()

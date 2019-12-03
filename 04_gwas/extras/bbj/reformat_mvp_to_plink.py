from __future__ import print_function
from __future__ import division 
import scipy.stats as st
import os, sys
import gzip 
import numpy as np
#import rpy2
#from rpy2.robjects import numpy2ri
#numpy2ri.activate()
#from rpy2.robjects import Environment

# Create an R environment
#env = Environment()

# Bind in R the R vector to the symbol "x" and
# in that environment
#env['x'] = rnorm(100)

# Build a tuple of pairs (<argument name>, <argument>).
# Note that the argument is a symbol. R will resolve what
# object is associated to that symbol when the function
# is executed.

def rev_comp(allele):
    comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a'}
    return(''.join([comp[x] if x in comp else x for x in allele[::-1]]))

fin = open(sys.argv[1],'r').readlines()
bq = sys.argv[2]
if bq == 'binqt':
    fout = gzip.open(os.path.join(os.path.dirname(sys.argv[1]), 'bin', os.path.basename(sys.argv[1]) +  '.plink.gz'), 'wt')
else:
    fout = gzip.open(os.path.join(os.path.dirname(sys.argv[1]), bq, os.path.basename(sys.argv[1]) +  '.plink.gz'), 'wt')

for l in fin[0:]:
    line = l.rstrip().split('\t')
    if len(line) < 22:
        continue
    if line[3] == "P-value":
        if bq == "bin":
            #CHROM  POS     ID      REF     ALT     A1      FIRTH?  TEST    OBS_CT  OR      SE      Z_STAT  P
            print('#CHROM\tPOS\tID\tREF\tALT\tA1\tFIRTH?\tTEST\tOBS_CT\tOR\tSE\tZ_STAT\tP', file = fout)
        if bq == "binqt":
            #CHROM  POS     ID      REF     ALT     A1      FIRTH?  TEST    OBS_CT  OR      SE      Z_STAT  P
            print('#CHROM\tPOS\tID\tREF\tALT\tA1\tFIRTH?\tTEST\tOBS_CT\tOR\tSE\tZ_STAT\tP', file = fout)
        elif bq == "qt":
            #CHROM  POS     ID      REF     ALT     A1      TEST    OBS_CT  BETA    SE      T_STAT  P
            print('#CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP', file = fout)
    else:
        base_index = 1 # 1 for 1-based index (in most cases), 0 for 0-based index
        if( base_index != 1 ):
            line[7] = str(int(line[7]) + 1 - base_index)  
        if( line[10] == '-' ):            # reverse complement
            for colidx in [11, 12, 19]:
                line[colidx] = rev_comp(line[colidx])

        if len(line[17]) < 1:
            if bq == "binqt":
                line[17] = float(line[16])/(-st.norm.ppf(float(line[3])))
            else:
                line[17] = np.log(float(line[16]))/(-st.norm.ppf(float(line[3])))
            if bq == "binqt":
                print(line[6],line[7],line[2], line[11], line[12], line[19],"N",  "ADD", line[20], float(line[16]), line[17], np.log(float(line[16]))/float(line[17]), line[3], sep = '\t',  file = fout)
            elif bq == "bin":
                print(line[6],line[7],line[2], line[11], line[12], line[19],"N",  "ADD", line[20], np.exp(float(line[16])), line[17], float(line[16])/float(line[17]), line[3], sep = '\t',  file = fout)
            elif bq == "qt":
                print(line[6],line[7],line[2], line[11], line[12], line[19], "ADD", line[20], line[16], line[17], float(line[16])/float(line[17]), line[3], sep = '\t',  file = fout)
        else:
            if bq == "binqt":
                print(line[6],line[7],line[2], line[11], line[12], line[19],"N",  "ADD", line[20], float(line[16]), line[17], np.log(float(line[16]))/float(line[17]), line[3], sep = '\t',  file = fout)
            elif bq == "bin":
                print(line[6],line[7],line[2], line[11], line[12], line[19],"N",  "ADD", line[20], np.exp(float(line[16])), line[17], float(line[16])/float(line[17]), line[3], sep = '\t',  file = fout)
            elif bq == "qt":
                print(line[6],line[7],line[2], line[11], line[12], line[19], "ADD", line[20], line[16], line[17], float(line[16])/float(line[17]), line[3], sep = '\t',  file = fout)

fout.close()


from __future__ import print_function

_README_ = '''
-------------------------------------------------------------------------
Check and make sure the REF/ALT alleles in GWAS results corresponds to A1/A2 alleles in UKBB bim file

Author: Yosuke Tanigawa (ytanigaw@stanford.edu)
Date: 2017/08/19, revised on 2017/11/12, 2019/3/29
-------------------------------------------------------------------------
'''
import sys, os, argparse, fileinput


import numpy as np
import pandas as pd


def flipfix_A1A2(in_f):
    in_file = '/dev/stdin' if in_f is None else in_f

    for line in fileinput.input(in_file):
        line_s = line.rstrip()
        l = line_s.split('\t')
        if(line_s[0] == '#'):
            col_idx_dict = dict(zip(l, np.arange(len(l))))
            # print header line
            print(line_s)
        else:
            CHR=str(l[col_idx_dict['#CHROM']])
            POS=str(l[col_idx_dict['POS']])
            REF=l[col_idx_dict['REF']]
            ALT=l[col_idx_dict['ALT']]
            A1=l[col_idx_dict['A1']]
            if(A1 == ALT):
                # if there is no flip
                print(line_s)
            elif(A1 == REF):
                # if there is a flip
                l[col_idx_dict['REF']] = ALT
                l[col_idx_dict['ALT']] = REF
                if('BETA' in col_idx_dict and l[col_idx_dict['TEST']] == 'ADD' and l[col_idx_dict['BETA']] != 'NA'):
                    l[col_idx_dict['BETA']] = -float(l[col_idx_dict['BETA']])
                elif('OR' in col_idx_dict and l[col_idx_dict['TEST']] == 'ADD' and l[col_idx_dict['OR']] != 'NA'):
                    l[col_idx_dict['OR']] = np.exp(- np.log(float(l[col_idx_dict['OR']])))
                print('\t'.join([str(x) for x in l]))
            else:
                # multi-allelic, haven't implemented yet for this case. Dump error and exit
                raise RuntimeError('Unexpected entry: {}'.format(line.rstrip()))
                sys.exit(1)

def main():  
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=_README_
    )    
    parser.add_argument('-i', metavar='i', default=None, help='input file [default: stdin]')
    args = parser.parse_args()
    flipfix_A1A2(args.i)
   
if __name__ == "__main__":
    main()


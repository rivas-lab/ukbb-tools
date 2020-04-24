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

def flipfix_A1A2(in_f, variants_file, to38):
    in_file = '/dev/stdin' if in_f is None else in_f

    for line in fileinput.input(in_file):
        line_s = line.rstrip()
        l = line_s.split('\t')
        if(line_s[0] == '#'):
            col_idx_dict = dict(zip(l, np.arange(len(l))))
            # print header line
            print(line_s)
            if to38 or 'A1' not in col_idx_dict:
                var_df = pd.read_csv(variants_file, sep='\t')
                dictA1 = dict(zip(var_df['ID'], var_df['hg37_A1']))
                dictA2 = dict(zip(var_df['ID'], var_df['hg37_A2']))
                dict38 = dict(zip(var_df['ID'], var_df['hg38_pos']))
        else:
            #CHR=str(l[col_idx_dict['#CHROM']])
            #POS=str(l[col_idx_dict['POS']])
            ID=l[col_idx_dict['ID']]
            REF=l[col_idx_dict['REF']]
            ALT=l[col_idx_dict['ALT']]
            if 'A1' in col_idx_dict:
                A1=l[col_idx_dict['A1']]
                A2=None
            else:
                # output from old version of PRISM
                A1=dictA1[ID]
                A2=dictA2[ID]
            if to38:
                # change coordinates
                l[col_idx_dict['POS']] = dict38[ID]

            if(  (A1 == ALT) and ((A2 is None) or (A2 == REF))):
                # if there is no flip
                print('\t'.join([str(x) for x in l]))
            elif((A1 == REF) and ((A2 is None) or (A2 == ALT))):
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
    parser.add_argument(
        '--variants', metavar='v', default='/oak/stanford/groups/mrivas/private_data/ukbb/24983/cal/pgen/ukb24983_cal_cALL_v2.liftOver.tsv.gz', 
        help='liftOver file'
    )
    parser.add_argument('--to38', action='store_true', help='change coordinates to hg38')
    args = parser.parse_args()
    flipfix_A1A2(args.i, args.variants, args.to38)
   
if __name__ == "__main__":
    main()


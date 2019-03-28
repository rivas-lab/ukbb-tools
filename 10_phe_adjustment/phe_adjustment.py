from __future__ import print_function
_README_='''
Apply phenotype adjustment
'''
import sys, os, argparse, fileinput
import numpy as np
import pandas as pd

def phe_adjustment(in_f, adjustment_factor, adjustment_phe):
    statin = pd.read_csv(
        adjustment_phe, sep='\t', usecols=[1,2]
    )
    statin_users = set([str(x) for x in statin[statin['statin'] == 2]['IID']])

    in_file = '/dev/stdin' if in_f is None else in_f
    for line in fileinput.input(in_file):
        line_s = line.rstrip()
        l = line_s.split('\t')
        if str(l[1]) not in statin_users:
            print('\t'.join([str(x) for x in l]))
        else:
            l[2] = float(l[2]) / adjustment_factor
            print('\t'.join([str(x) for x in l]))

def main():  
    _default_adjustment_phe=os.path.join(
        os.path.dirname(os.path.realpath(__file__)), 
	'statin_phes', 'statin.phe'
    )
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=_README_
    )    
    parser.add_argument('-i', metavar='i', default=None, help='input file [default: stdin]')
    parser.add_argument('--mult', metavar='m', type=float, help='adjustment multiplier (apply * 1/m for cases)')
    parser.add_argument(
        '--adjustment_phe', metavar='p', default=_default_adjustment_phe,
        help='adjustment phe file'
    )
    args = parser.parse_args()
    phe_adjustment(in_f = args.i, adjustment_factor = args.mult, adjustment_phe = args.adjustment_phe)
   
if __name__ == "__main__":
    main()


from __future__ import print_function
_README_ = '''
-------------------------------------------------------------------------
Prepare phe files from provided UKB Field/Coding specs.

For binary traits: specify the coding of cases and controls. Others will be treated as missing.

For quantitative traits: (optionally) specify ordering of phenotype values, or values to exclude.
 
Authors: Yosuke Tanigawa (ytanigaw@stanford.edu) and Matthew Aguirre (magu@stanford.edu)
Date: 2018/02/13
-------------------------------------------------------------------------
'''

import numpy as np
import pandas as pd
import collections as cl
import argparse
import datetime
from functools import reduce
import os, sys

# logging

import logging
from logging.config import dictConfig
dictConfig(dict(
    version = 1,
    formatters = {
        'f': {'format':
              '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'}
        },
    handlers = {
        'h': {'class': 'logging.StreamHandler',
              'formatter': 'f',
              'level': logging.DEBUG}
        },
    root = {
        'handlers': ['h'],
        'level': logging.DEBUG,
        },
))


def get_tsv_from_tab(in_tab, field_id):
    return pd.read_csv(in_tab, sep='\t', usecols=lambda col: col=="f.eid" or field_id == col.split(".")[1])


def create_bin_phe_file(in_tsv, out_phe, out_log, field_id, case, control, missing_is_control):
    logger_phe = logging.getLogger('create_phe_file')    
    hdlr = logging.FileHandler(out_log)
    formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    hdlr.setFormatter(formatter)    
    logger_phe.addHandler(hdlr) 
    
    tsv_df=get_tsv_from_tab(in_tsv, field_id)
    logger_phe.info('encoding: case: ' + ';'.join([str(x) for x in sorted(case)]))
    logger_phe.info('encoding: control: ' + ';'.join([str(x) for x in sorted(control)])    )
    logger_phe.info('number_of_measurements: ' + str(tsv_df.shape[1] - 1))
    logger_phe.info('measurement_names: ' + str(';'.join(tsv_df.columns[1:])))

    # extract data 
    data = np.array(tsv_df.iloc[:,1:], dtype=np.float64)
   
    # convert coding into 
    data_case    = np.vectorize(lambda x: (not np.isnan(x)) and (str(int(x)) in case))(data)    
    data_control = np.vectorize(lambda x: (not np.isnan(x)) and (str(int(x)) in control))(data)    
    data_missing = np.invert( data_case | data_control )

    logger_phe.info('number_of_cases: ' + ';'.join([str(x) for x in np.sum(data_case, axis=0)]))
    logger_phe.info('number_of_controls: ' + ';'.join([str(x) for x in np.sum(data_control, axis=0)]))
    logger_phe.info('number_of_missing_values: ' + ';'.join([str(x) for x in np.sum(data_missing, axis=0)])) 

    # aggregate multiple measurments
    data_agg_case    = np.apply_along_axis(lambda l: reduce(lambda x, y: x or y, l), 1, data_case)
    data_agg_missing = np.apply_along_axis(lambda l: reduce(lambda x, y: x and y, l), 1, data_missing)
    data_agg_control = np.invert( data_agg_case | data_agg_missing )

    logger_phe.info('number_of_cases_(aggregated): {}'.format(np.sum(data_agg_case, axis=0)))
    logger_phe.info('number_of_controls_(aggregated): {}'.format(np.sum(data_agg_control, axis=0)))
    logger_phe.info('number_of_missing_values_(aggregated): {}'.format(np.sum(data_agg_missing, axis=0))) 
 
    data_agg = (np.where(data_agg_case, 2, 0) + np.where(data_agg_control, 1, 0) + np.where(data_agg_missing, -9, 0))
    # re-encode missing values, if applicable
    if missing_is_control:
        np.place(data_agg, data_agg == -9, 1)
    
    # write to a .phe file
    phe_df = pd.DataFrame({
        'FID'  : tsv_df.ix[:, 0],
        'IID'  : tsv_df.ix[:, 0],
        'data' : data_agg
    }).sort_values('FID').to_csv(out_phe, sep='\t', index=False, header=False)    
    return


def create_qt_phe_file(in_tsv, out_phe, out_log, field_id, order=[], exclude=[]):
    logger_phe = logging.getLogger('create_phe_file')    
    hdlr = logging.FileHandler(out_log)
    formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    hdlr.setFormatter(formatter)    
    logger_phe.addHandler(hdlr) 
    
    # extract data  
    tsv_df=get_tsv_from_tab(in_tsv, field_id)
    data_raw = np.array(tsv_df.iloc[:,1:], dtype=np.float64)
    
    logger_phe.info('number_of_measurements: ' + str(tsv_df.shape[1]-1))
    logger_phe.info('measurement_names: ' + str(';'.join(tsv_df.columns[1:])))
    
    if(exclude[0] != ""):
        exclude=sorted([float(u) for u in exclude])
        logger_phe.info('excluding these values: ' + str(';'.join([str(x) for x in exclude])))
        data_raw = np.where(np.any(np.isnan(np.array([np.where(data_raw==i,np.nan,data_raw) for i in exclude])),axis=0),np.nan,data_raw)
    
    if(order[0] != ""): # remap values, if applicable
        order=sorted([float(u) for u in order]) # some tables are split by commas instead ):
        logger_phe.info('order of coded values: ' + str('<'.join([str(int(x)) for x in order])))
        data_decode_dict = dict(zip(order,range(1, len(order) + 1)))
        data = np.array([np.nan if np.isnan(x) else data_decode_dict[x] if x in data_decode_dict else np.nan for row in data_raw for x in row]).reshape(data_raw.shape)
    else:
        data = data_raw
    print(data.shape) 
    logger_phe.info('number_of_nans: ' + ';'.join([str(x) for x in np.sum(np.isnan(data), axis=0)]))

    # take median
    aggregated = np.nanmedian(data, axis=1)

    logger_phe.info('final_number_of_nan: ' + str(np.sum(np.isnan(aggregated))))
    
    # write to a .phe file
    phe_df = pd.DataFrame({
        'FID' : tsv_df.ix[:, 0],
        'IID' : tsv_df.ix[:, 0],
        'data': ['-9' if np.isnan(x) else str(x) for x in aggregated]
    }).to_csv(out_phe, sep='\t', index=False, header=False)
    return
 

def main():    
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=_README_
    )    
    parser.add_argument('-i', metavar='i', required=True,
                        help='input tab file')        
    parser.add_argument('-f', metavar='f', required=True,
                        help='input field id')        
    parser.add_argument('-o', metavar='o', required=True,
                        help='output phe file')
    parser.add_argument('-l', metavar='l', required=True,
                        help='output log file')
    parser.add_argument('--type', metavar='bin|qt', dest='t', required=True,
                        help='phe type (must be bin|qt)')
    parser.add_argument('--case', metavar='2', type=str, nargs='*',
                        help='comma delimited coding for cases')
    parser.add_argument('--control', metavar='C', type=str, nargs='*',
                        help='comma delimited coding for controls')
    parser.add_argument('--order', metavar='O', dest='O', nargs='*', type=str,
                        help='semicolon delimited coding order (increasing) for QTs')
    parser.add_argument('--exclude', metavar='X', dest='X', nargs='*', type=str,
                        help='semicolon delimited values to exclude for QTs -- must be sanitized by being surrounded in "-quotes"')
    parser.add_argument('--missing-controls', dest='M', action='store_true', 
                        help='flag to expand control set to include missing values as controls')
     
    args = parser.parse_args()
    print(args)
    if args.t == 'bin':
        if (not args.case) or (not args.control and not args.M): 
            sys.stderr.write('['+os.path.abspath(sys.argv[0]) + 
                             datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S") + 
                             '] --case and --control are required for binary phenotypes!')
            sys.exit(1)
        create_bin_phe_file(
            in_tsv   = args.i, 
            out_phe  = args.o, 
            out_log  = args.l, 
            field_id = args.f,
            # splitting and replacing like this enable: (1) ";" and "," delimiters, and 
            #                                           (2) inputs sanitized with quotes (in case of leading minus sign) 
            case     = set(';'.join(args.case).replace('"','').replace("'","").replace(',',';').split(';')), 
            control  = set(';'.join(args.control).replace('"','').replace("'","").replace(',',';').split(';')),
            missing_is_control  = args.M 
        )
    elif args.t == 'qt':
        create_qt_phe_file(
            in_tsv   = args.i, 
            out_phe  = args.o, 
            out_log  = args.l, 
            field_id = args.f,
            # these (and case/control above) return [''] if args.X is empty or is ['']
            exclude  = ';'.join(args.X).replace('"','').replace("'","").replace(',',';').split(';') if args.X else [''], 
            order    = ';'.join(args.O).replace('"','').replace("'","").replace(',',';').split(';') if args.O else ['']
        )
    else:
        sys.stderr.write('['+os.path.abspath(sys.argv[0]) + 
                         datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S") + 
                         '] unsupported phe type: '+args.t)
        sys.exit(1)
    
    sys.stderr.write('['+os.path.abspath(sys.argv[0]) +
                     datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S") +
                     '] phe file is written to: '+args.o + '\n')
    return

if __name__ == "__main__":
    main()


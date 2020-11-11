import sys


import pandas as pd
import pyarrow
from pyarrow import csv, parquet, feather


def get_dtype_master_gwas():
    '''
    Here, we define the data types in the master GWAS file
    '''
    return {
        '#CHROM'           : 'str',
        'CHROM'            : 'str',
        'POS'              : 'int64',
        'Variant_ID'       : 'str',
        'Csq'              : 'str',
        'GBE_ID'           : 'str',
        'population'       : 'str',
        'REF'              : 'str',
        'ALT'              : 'str',
        'A1'               : 'str',
        'TEST'             : 'str',
        'OBS_CT'           : 'int64',
        'BETA'             : 'float64',
        'SE'               : 'float64',
        'T_or_Z_STAT'      : 'float64',
        'P'                : 'float64',
        'Direction'        : 'str',
        'HetISq'           : 'float64',
        'HetChiSq'         : 'float64',
        'HetDf'            : 'int32',
        'HetPVal'          : 'float64',
        'ERRCODE'          : 'str',
        'FIRTH?'           : 'str',
        'OR'               : 'float64',
        'SYMBOL'           : 'str',
        'GBE_short_name'   : 'str',
        'Consequence'      : 'str',
        'Gene'             : 'str',
        'GBE_category'     : 'str'
    }


def csv_to_parquet(in_f, out_f, delimiter='\t', dtype=None):
    '''
    Read a csv file as a stream and save it as Apache Parquet file
    '''
    pa_reader = pyarrow.csv.open_csv(
        in_f,
        # read_options = pyarrow.csv.ReadOptions(use_threads=True),
        read_options = pyarrow.csv.ReadOptions(use_threads=True, skip_rows=1, column_names=[x.replace('#', '') for x in pd.read_csv(in_f, sep='\t', nrows=0).columns]),
        parse_options=pyarrow.csv.ParseOptions(delimiter=delimiter),
        convert_options=pyarrow.csv.ConvertOptions(column_types=dtype)
    )
    # when the header line in the input line starts with '#', 
    # we can specify the following read_options:
    #   pyarrow.csv.ReadOptions(
    #       use_threads=True,
    #       skip_rows=1, # here we provide `skip_rows` and `column_names` to get rid of '#'
    #       column_names=[x.replace('#', '') for x in pd.read_csv(in_f, sep='\t', nrows=0).columns]
    #   )
    
    pa_writer = pyarrow.parquet.ParquetWriter(
        out_f, pa_reader.schema, compression='zstd'
    )

    nrow=0
    for batch in pa_reader:
        batch_df=batch.to_pandas()
        nrow += batch_df.shape[0]
        pa_writer.write_table(pyarrow.Table.from_pandas(batch_df))
        # in principle, it should be possible to directly save it without converting to/from pandas df
        # however, when I try 
        #   pa_writer.write_table(pyarrow.Table.from_batches(batch))
        # I got and error:
        #   TypeError: Cannot convert pyarrow.lib.StringArray to pyarrow.lib.RecordBatch
    return nrow

def csv_to_feather(in_f, out_f, delimiter='\t', dtype=None):
    '''
    Read a csv file as a stream and save it as Apache Feather file
    It seems like the stream writer has not implemented yet
    '''
    pa_reader = pyarrow.csv.open_csv(
        in_f,
        read_options = pyarrow.csv.ReadOptions(use_threads=True),
        parse_options=pyarrow.csv.ParseOptions(delimiter=delimiter),
        convert_options=pyarrow.csv.ConvertOptions(column_types=dtype)
    )
    df = pa_reader.read_pandas()
#     df = pd.read_csv(
#         in_f, sep=delimiter, dtype=dtype, engine='c'
# #     ).rename(
# #         columns = {'#CHROM' : 'CHROM'}
#     )    
    pyarrow.feather.write_feather(
        df,
        out_f,
        compression='zstd'
    )
    return df.shape[0]

def tsv2columnar(in_f, out_f, columnar, delimiter='\t', dtype=None):
    if(columnar == 'parquet'):
        csv_to_parquet(in_f, out_f, delimiter=delimiter, dtype=dtype)
    elif(columnar == 'feather'):
        csv_to_feather(in_f, out_f, delimiter=delimiter, dtype=dtype)

if __name__ == '__main__':
    in_f     = sys.argv[1]
    out_f    = sys.argv[2]
    columnar = sys.argv[3]
    tsv2columnar(in_f, out_f, columnar, delimiter='\t', dtype=get_dtype_master_gwas())

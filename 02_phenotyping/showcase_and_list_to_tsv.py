_README='read in list of fields, subset from big csv.'


import argparse


import pandas as pd


def join_and_add_cols(filename, ref):
    dds = pd.read_csv(ref)

    #Read in the fields
    fields = pd.read_csv(filename, header=None, names=['FieldID'])

    #Subset the df with a merge
    subsetted = fields.merge(dds, left_on='FieldID', right_on='FieldID')

    #Add in the extra column names
    add_colnames = ['Annotator', 'Annotation date', 'Name', 'GBE ID', 'TableID', 'Field', 'QT_total_num', 'BIN_total_num', 'QT_index', 'BIN_index', 'coding_exclude', 'coding_QT', 'coding_binary_case', 'coding_binary_control']
    subsetted = pd.concat([subsetted,pd.DataFrame(columns=add_colnames)], sort=True)

    #Reorder the columns
    new_col_order = ['Annotator', 'Annotation date', 'Name', 'GBE ID', 'TableID', 'Field', 'FieldID', 'QT_total_num', 'BIN_total_num', 'QT_index', 'BIN_index', 'coding_exclude', 'coding_QT', 'coding_binary_case', 'coding_binary_control', 'Participants', 'Stability', 'ValueType', 'Units', 'Strata', 'Sexed', 'Instances', 'Array', 'Coding', 'Link']
    subsetted = subsetted[new_col_order]

    #ensure types of columns are int
    subsetted = subsetted.astype({"FieldID": int, "Participants": int, "Instances": int, "Array": int})
    
    return(subsetted)

  
def showcase_and_list_to_tsv_main():
    #Parse args to get the csv of fields
    parser = argparse.ArgumentParser(description=_README)
    parser.add_argument(
      'csvfile', 
      type=argparse.FileType('r'), 
      help='Input csv file'
    )
    parser.add_argument(
      '--ref', type=argparse.FileType('r'), 
      help='Data_Dictionary_Showcase.csv',
      default='Data_Dictionary_Showcase.csv'
    )
    parser.add_argument(
      '--out', 
      help='output tsv file (default: ${csvfile%.csv}.tsv)',
      default=None
    )
    
    args = parser.parse_args()

    filename = args.csvfile
    ref = args.ref
    if args.out is not None:
        outfilename = args.out
    else:
        outfilename = filename.name[:-3] + 'tsv'

    subsetted = join_and_add_cols(filename, ref)
    
    print(outfilename)
    #Export to new file
    subsetted.to_csv(outfilename, sep='\t', index=False)

    
if __name__ == '__main__':
    showcase_and_list_to_tsv_main()

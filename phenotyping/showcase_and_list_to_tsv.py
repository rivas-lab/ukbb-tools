import pandas as pd
import argparse

dds = pd.read_csv('Data_Dictionary_Showcase.csv')

#Parse args to get the csv of fields
parser = argparse.ArgumentParser(description='read in list of fields, subset from big csv.')
parser.add_argument('csvfile', type=argparse.FileType('r'), help='Input csv file')
args = parser.parse_args()

#Read in the fields
fields = pd.read_csv(args.csvfile, header=None, names=['FieldID'])

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

#Export to new file
subsetted.to_csv('google_sheet_upload.tsv', sep='\t', index=False)

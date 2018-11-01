#!/bin/python
import os

__README__="""


Author: Matthew Aguirre (SUNET: magu)
"""

in_tsv="/oak/stanford/groups/mrivas/dev-ukbb-tools/phenotypes/tables/ukb9797_20171211.tsv"
ukb_table_id="9797" # this can be in the table as well -- could be arbitrarily specified 
out_dir="/oak/stanford/groups/mrivas/dev-ukbb-tools/phenotypes/9797"

# the below columns are zero-indexed in the input tsv
name_col=0  # column for output name (gbe id)
field_id=2  # column for ukb field id
case_col=29 # column for binary case values
ctrl_col=30 # column for binary control values
qtfc_col=28 # column for ordering of quantitative values, least to greatest (from categorical field)
excl_col=27 # column for values to mark as missing
# todo: missing as control?

whole_table=True # define whole table or split out into new jobs?
header=True

if whole_table:
	from make_phe import *
	phe_info = {}
	with open(in_tsv, 'r') as f:
		for i,line in enumerate(f):
			if header and i==0:
				continue
			fields = line.rstrip().split('\t')
			if not fields[name_col]: continue
			# the control switch is in case empty cells get pruned from the TSV
			phe_info[ fields[name_col] ] = {'case':     fields[case_col] if case_col < len(fields) else '',
							'control':  fields[ctrl_col] if ctrl_col < len(fields) else '',
							'qt_order': fields[qtfc_col] if qtfc_col < len(fields) else '',
							'exclude':  fields[excl_col] if excl_col < len(fields) else '',
							'field_id': fields[field_id] if field_id < len(fields) else ''}
	# this info should be logged somewhere
	for phe_name, phe_values in phe_info.items():
		print(phe_name, phe_values)
		if phe_values['case']: # assume binary if we have a case definition
			create_bin_phe_file(in_tsv=os.path.join("/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/download",
								   str(int(ukb_table_id)-1), 
								   ukb_table_id, 
								   "ukb{}.tab".format(ukb_table_id)),
					       out_phe=os.path.join(out_dir, "{0}.phe".format(phe_name)),
					       out_log=os.path.join(out_dir, "logs", "{0}.log".format(phe_name)),
					       field_id=phe_values['field_id'],
					       case=phe_values['case'].replace(',',';').split(';'),
					       control=phe_values['control'].replace(',',';').split(';'),
					       missing_is_control=False) # todo)
		else: # assume qt
			create_qt_phe_file(in_tsv=os.path.join("/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/download",
							       str(int(ukb_table_id)-1), 
							       ukb_table_id, 
							       "ukb{}.tab".format(ukb_table_id)),
					   out_phe=os.path.join(out_dir, "{0}.phe".format(phe_name)),
					   out_log=os.path.join(out_dir, "logs", "{0}.log".format(phe_name)),
					   field_id=phe_values['field_id'],
					   order=phe_values['qt_order'].replace(',',';').split(';'),
					   exclude=phe_values['exclude'].replace(',',';').split(';'))
else:
	pass
	# todo: split out into individual jobs (presumably to run gwas, so not sure i want to keep this as an option

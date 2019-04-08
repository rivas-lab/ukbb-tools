#!/bin/python
from math import exp,log

with open('/oak/stanford/groups/mrivas/ukbb/master_phe/decomposition_phe_20171220.phe', 'r') as f, open('INI1003064.phe', 'w') as o:
	for i,line in enumerate(f):
		x = line.split()
		if not i:
			a,b,c,d = (x.index('age'), x.index('sex'), x.index('INI50'), x.index('INI3064'))
			print a,b,c,d
		else:
			# catch error for bad definition in master
			try:
				age, sex, height, pef = float(x[a]), float(x[b]) + 1, float(x[c]), float(x[d])
			except IndexError:
				o.write('\t'.join([x[0], x[0], '-9']) + '\n')
				continue
			# catch error for missing/bad values in master
			if (not all([age, sex, height, pef])) or any([height < 0, height=='NA', pef < 0, pef=='NA']):
				o.write('\t'.join([x[0], x[0], '-9']) + '\n')
				continue
			# case management for sex, and catch misencoded sex
			if sex == 1:
				o.write('\t'.join([x[0], x[0], str(pef/exp(0.544*log(float(age)) - 0.0151*float(age) - 74.7/height + 5.48))]) + '\n')
			elif sex == 2:
				o.write('\t'.join([x[0], x[0], str(pef/exp(0.376*log(float(age)) - 0.012*float(age) - 58.8/height + 5.63))]) + '\n')
			else:
				o.write('\t'.join([x[0], x[0], '-9']) + '\n')


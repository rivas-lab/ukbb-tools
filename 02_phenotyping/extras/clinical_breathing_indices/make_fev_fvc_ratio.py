#!/bin/python
from math import exp,log

with open('/oak/stanford/groups/mrivas/ukbb/master_phe/decomposition_phe_20171220.phe', 'r') as f, open('INI1003063.phe', 'w') as o:
	for i,line in enumerate(f):
		x = line.split()
		if not i:
			a,b = x.index('INI3063'), x.index('INI3062')
		else:
			# catch error in master phe
			try:
				fev, fvc = float(x[a]), float(x[b])
			except IndexError:
				o.write('\t'.join([x[0], x[0], '-9']) + '\n')
				continue
			# catch bad/missing phenotype value
			if (not all([fev, fvc])) or any([fev < 0, fvc < 0]):
				o.write('\t'.join([x[0], x[0], '-9']) + '\n')
				continue
			o.write('\t'.join([x[0], x[0], str(fev/fvc)]) + '\n')

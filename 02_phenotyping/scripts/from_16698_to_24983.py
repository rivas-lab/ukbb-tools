import sys

# usage: python to24983.py in_file out_file
in_file=sys.argv[1]
out_file=sys.argv[2]

# load ref, key on 16698:
with open('/oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_16698_mapping.tsv', 'r') as f:
    # lines are 24983_id 16698_id
    to24983 = dict([line.split()[::-1] for line in f])

# walk through in_file, replace id and write to out_file
with open(in_file, 'r') as i, open(out_file, 'w') as o:
    for line in i:
        split_line = line.split()
        if split_line[0] in to24983:
            o.write('\t'.join([to24983[item] if n < 2 else item for n, item in enumerate(split_line)])+'\n')

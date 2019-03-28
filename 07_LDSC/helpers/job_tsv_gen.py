import itertools
with open('sumstats.lst') as f:
    l = f.read().splitlines()
with open('jobs.tsv', 'w') as f:
    for x, y in itertools.combinations(l, 2):
        f.write('{}\t{}\n'.format(x, y))


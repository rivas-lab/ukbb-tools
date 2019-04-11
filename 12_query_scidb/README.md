# Querying SciDB server to get the PheWAS results

To generate PheWAS plot using the data on GBE/SciDB, one can use `phewas_query.py` to export the data in SciDB.

Usage of the script is as follows.

```
$ python /home/$USER/repos/rivas-lab/ukbb-tools/12_query_scidb/phewas_query.py -h
usage: phewas_query.py [-h] -i i -o o [--port P] [-n N] [-p p] [-l L]

-------------------------------------------------------------------------
phewas.py

Query SciDB (GBE) and get the PheWAS sumstats

Author: Yosuke Tanigawa (ytanigaw@stanford.edu)
Date: 2019/04/11
-------------------------------------------------------------------------

optional arguments:
  -h, --help  show this help message and exit
  -i i        query (comma-separated coordinates like 5-145895394,11-14865399)
  -o o        output npz file
  --port P    SciDB port (default: 8080)
  -n N        minimum case count (default: 100)
  -p p        p-value threshold (default: 0.001)
  -l L        a file that has a list of phenotypes (GBE_IDs)
```

In practical example, one may run this script like this:

```
$ python phewas_query.py -i 5-145895394,11-14865399 -o sumstats.out.1.tsv
```

If you are querying for multiple variants, it might be useful to get sumstats for the common set of phenotypes.
In that case, you can run the followings:

```
$ python phewas_query.py -i 5-145895394,11-14865399 -o sumstats.out.1.tsv
$ sumstats.out.1.tsv | awk '(NR>1){print $2}' > sumstats.out.1.phe.lst
$ python phewas_query.py -i 5-145895394,11-14865399 -l sumstats.out.1.phe.lst -p 1 -o sumstats.out.2.tsv
```

Note that this script depends on the SciDB process, SciDB-python package, and GBE source code (for schema of the SciDB table).


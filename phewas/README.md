## PheWAS

This directory contains scripts for running PheWAS (`phewas.py`) and generating the files necessary to do so (`combine_phe.py`).

# Dependencies:

_Note_: this section only needs to be updated when new phenotype data comes in, but feel free to read it for more info about where the files used by `phewas.py` come from.

In order to run PheWAS, we need a master table containing all the phenotype data we've processed (see the [phenotyping](https://github.com/rivas-lab/ukbb-tools/blob/master/phenotyping/) folder for more info). To do this, we leverage the phenotype info table, which is also in the phenotyping folder. The maker script, `combine_phe.py` will walk through this file and selecting phenotypes with n/N > 100 for inclusion. It then loads each of those files into memory (!!) via a dictionary keyed on sample ID, then writes all the data out to file. 

If you're rerunning this, you'll probably want a lot of memory -- 64g is suffient in my experience (Matt).


# Analysis:

This script is under construction. 

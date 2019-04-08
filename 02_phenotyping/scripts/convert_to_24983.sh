#/bin/bash
# replace dir in ls below with desired directory to convert to 24983 from 16698
# replace python script path if not proper
for f in $(ls /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/rohit/phe); do
    python /oak/stanford/groups/mrivas/users/$USER/repos/ukbb-tools/02_phenotyping/scripts/from_16698_to_24983.py $f to24983$f
    mv to24983$f $f
done

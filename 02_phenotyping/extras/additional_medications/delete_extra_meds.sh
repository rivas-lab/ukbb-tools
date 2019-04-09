#/bin/bash

names=/oak/stanford/groups/mrivas/users/guhan/repos/ukbb-tools/02_phenotyping/extras/additional_medications/med_names.txt

for f in $(ls | grep -v delete); do
    base=$(basename "$f" .phe)
    if ! grep -q $base "$names"
    then
        rm -rf $base.phe
    fi
done

#!/bin/bash
set -beEuo pipefail

GBE_IDs=(
MED1140861958
MED1140888594
MED1140888648
MED1141146234
MED1141192410
MED1141146138
)

# this one was missing from the master phe file
# https://github.com/rivas-lab/ukbb-tools/issues/7
echo MED1140861922
if [ ! -f  statin_phes/MED1140861922.phe ] || [ -z  statin_phes/MED1140861922.phe ] ; then
{
echo "#FID IID MED1140861922"
bash /oak/stanford/groups/mrivas/ukbb24983/ukb24983_16698_mapping.sh \
/oak/stanford/groups/mrivas/private_data/ukbb/16698/phe/ukb9797_20170827_20003/MED1140861922.phe 
} | tr " " "\t" > statin_phes/MED1140861922.phe
fi

for i in ${GBE_IDs[@]} ; do
    echo $i
    if [ ! -f  statin_phes/$i.phe ] || [ -z  statin_phes/$i.phe ] ; then
        bash ../05_phewas/extract_phe.sh $i > statin_phes/$i.phe
    fi
done

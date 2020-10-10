#!/bin/bash
set -beEuo pipefail

unfinished=2c_find_unfinished.$(date +%Y%m%d-%H%M%S).lst
oom=${unfinished%.lst}.oom.lst
others=${unfinished%.lst}.others.lst

out_d=/scratch/groups/mrivas/ukbb24983/cnv/annotation_20201003/output_vep

find ${out_d} -name "*.tsv" | awk -v FS='.' '{print $5}' | sed -e 's/^00//g' | sed -e 's/^0//g' | sort | comm -3 <(seq 998 | sort) /dev/stdin | sort -n > ${unfinished}

grep -l 'oom-kill event' logs/vep.9222768_*.err  logs/vep.9222540_*.err \
    | awk -v FS='_' '{print $NF}' | awk -v FS='.' '{print $1}' | sort | comm -12 <(cat ${unfinished} | sort) /dev/stdin | sort -n > ${oom}

cat ${unfinished} | sort | comm -3 /dev/stdin <(cat ${oom} | sort ) | sort -n > ${others}


echo "oom"
cat ${oom} | tr '\n' ',' | rev | cut -c2- | rev

echo "others"
cat ${others} | tr '\n' ',' | rev | cut -c2- | rev

echo "clean-up"
cat ${oom} \
    | while read idx ; do 
    idx_pad=$(perl -e "print(sprintf('%03d', ${idx}))")
    find ${out_d} -name "*\.${idx_pad}\.*" | sort 
done | while read f ; do 
    mv $f $(dirname ${out_d})/output_vep_bkup/
done


#!/bin/bash
set -beEuo pipefail

cat 2d_clean_zstd_conflicts.20201104-142033.txt | egrep -v '.log$' | grep 'batch' | sed -e 's/.zst//g' | sort --parallel 6 | uniq -c \
    | awk '(NR>1){print $2}'

exit 0
############
bash 2e_list_zstd_duplicates.sh >  2e_list_zstd_duplicates.$(date +%Y%m%d-%H%M%S).txt


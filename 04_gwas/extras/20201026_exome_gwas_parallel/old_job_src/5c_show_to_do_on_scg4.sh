#!/bin/bash
set -beEuo pipefail
#finished_list=5b_finished_list.20201029-193622.tsv
finished_list=5b_finished_list.20201101-103415.tsv

scg_sh2_lst=old_job_src/squeue.20201029-001018.sh2-help.GBE_ID.txt

cat GBE_IDs_copy/FH.lst GBE_IDs_copy/cancer.lst GBE_IDs_copy/BIN.lst | egrep -v '#' | sort | comm -23 /dev/stdin <(cat ${finished_list} | awk -v p=white_british '($1==p){print $2}' | cat ${scg_sh2_lst} /dev/stdin |  sort ) | sort -V


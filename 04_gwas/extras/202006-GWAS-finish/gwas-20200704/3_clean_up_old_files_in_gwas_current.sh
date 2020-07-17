#!/bin/bash
set -beEuo pipefail

check_file="2_gwas.jobs.status.20200716-232605.tsv"

cat ${check_file} |
awk -v FS='\t' -v OFS='\t' -v T="TRUE" 'NR==1 || ($5 == T && $6 == T){print $7}' |
while read f ; do
    if [ -f $f ] ; then
        # echo $(zcat $f | wc -l) $f
        bkup_d="old_files_in_gwas_current/$(basename $(dirname $f))"
        if [ ! -d ${bkup_d} ] ; then mkdir -p ${bkup_d} ; fi
        mv $f ${bkup_d}/
    fi
done

cat ${check_file} |
awk -v FS='\t' -v OFS='\t' -v T="TRUE" 'NR==1 || ($5 == T && $6 == T){print $8}' |
while read f ; do
    if [ -f $f ] ; then
        # echo $(zcat $f | wc -l) $f
        bkup_d="old_files_in_gwas_current-symlinks/$(basename $(dirname $f))"
        if [ ! -d ${bkup_d} ] ; then mkdir -p ${bkup_d} ; fi
        mv $f ${bkup_d}/
    fi
done

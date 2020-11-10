#!/bin/bash
set -beEuo pipefail

plink_gz=$1
out_f=$(dirname $(dirname ${plink_gz}))/wc_l/$(basename $(dirname ${plink_gz}))/$(basename ${plink_gz%.gz}).wc_l.txt

if [ $# -gt 1 ] ; then out_f=$2 ; fi
if [ ! -d $(dirname ${out_f}) ] ; then mkdir -p $(dirname ${out_f}) ; fi

cat_or_zcat () {
    local file=$1
    if [ "${file%.gz}.gz" == "${file}" ] || [ "${file%.bgz}.bgz" == "${file}" ] ; then 
        zcat ${file} 
    elif [ "${file%.zst}.zst" == "${file}" ] ; then 
        zstdcat ${file}
    else
        cat ${file}
    fi
}

{
    echo "#file NF_wc_l wc_l non_NA_wc_l"     
    echo ${plink_gz} $(cat_or_zcat ${plink_gz} | awk -v FS='\t' '{print NF}' | sort | uniq -c | wc -l) $(cat_or_zcat ${plink_gz} | wc -l) $(cat_or_zcat ${plink_gz} | awk 'NR>1 && $(NF-1) != "NA"' | wc -l)
} | tr ' ' '\t' > ${out_f}

exit 0
##############
20013
7a_wc_l.sh

sbatch -p mrivas,normal --time=12:0:00 --mem=8000 --nodes=1 --cores=1 \
    --job-name=wc_l --output=logs/wc_l.sh.%A_%a.out --error=logs/wc_l.sh.%A_%a.err \
    --array=1-401 ${parallel_sbatch_no_err_check_sh} 7a_wc_l.sh 7b_wc_l_joblist.lst 50

452
7b_wc_l_joblist.20201104-114736.lst

sbatch -p mrivas,normal --time=12:0:00 --mem=8000 --nodes=1 --cores=1 \
    --job-name=wc_l --output=logs/wc_l.sh.%A_%a.out --error=logs/wc_l.sh.%A_%a.err \
    --array=1-91 ${parallel_sbatch_no_err_check_sh} 7a_wc_l.sh 7b_wc_l_joblist.20201104-114736.lst 5

315
7b_wc_l_joblist.20201104-131414.lst

sbatch -p mrivas,normal --time=12:0:00 --mem=8000 --nodes=1 --cores=1 \
    --job-name=wc_l --output=logs/wc_l.sh.%A_%a.out --error=logs/wc_l.sh.%A_%a.err \
    --array=1-105 ${parallel_sbatch_no_err_check_sh} 7a_wc_l.sh 7b_wc_l_joblist.20201104-131414.lst 3

385
7b_wc_l_joblist.20201104-145250.lst

sbatch -p mrivas,normal --time=12:0:00 --mem=8000 --nodes=1 --cores=1 \
    --job-name=wc_l --output=logs/wc_l.sh.%A_%a.out --error=logs/wc_l.sh.%A_%a.err \
    --array=1-129 ${parallel_sbatch_no_err_check_sh} 7a_wc_l.sh 7b_wc_l_joblist.20201104-145250.lst 3

324
7b_wc_l_joblist.20201104-161740.lst

sbatch -p mrivas,normal --time=2:0:00 --mem=8000 --nodes=1 --cores=1 \
    --job-name=wc_l --output=logs/wc_l.sh.%A_%a.out --error=logs/wc_l.sh.%A_%a.err \
    --array=1-108 ${parallel_sbatch_no_err_check_sh} 7a_wc_l.sh 7b_wc_l_joblist.20201104-161740.lst 3

323
7b_wc_l_joblist.20201104-170548.lst

sbatch -p mrivas,normal --time=2:0:00 --mem=8000 --nodes=1 --cores=1 \
    --job-name=wc_l --output=logs/wc_l.sh.%A_%a.out --error=logs/wc_l.sh.%A_%a.err \
    --array=1-108 ${parallel_sbatch_no_err_check_sh} 7a_wc_l.sh 7b_wc_l_joblist.20201104-170548.lst 3

17
7b_wc_l_joblist.20201104-205648.lst

sbatch -p mrivas,normal --time=1:0:00 --mem=8000 --nodes=1 --cores=1 \
    --job-name=wc_l --output=logs/wc_l.sh.%A_%a.out --error=logs/wc_l.sh.%A_%a.err \
    --array=1-17 ${parallel_sbatch_no_err_check_sh} 7a_wc_l.sh 7b_wc_l_joblist.20201104-205648.lst 1

#!/bin/bash
set -beEuo pipefail

mem=16000
cores=4

pop=$1 # related | others
batch_idx=$2
if [ $# -gt 2 ] ; then
    regression_type=$3 # linear
else
    regression_type="logistic"
fi
batch_size_bin=4
batch_size_qt=10

is_scg=$(echo scg $(hostname) | tr " " "\n" | grep scg | wc -l)
avx2_flag=$( cat /proc/cpuinfo  | grep flags | uniq  | awk -v FS=':' '{print $2}' | tr " " "\n" | cat /dev/stdin <(echo avx2) | grep -i avx2 | wc -l)
if [ "${is_scg}" -eq 2 ] ; then
    export MODULEPATH="$HOME/.modules:/labs/mrivas/.modules:$MODULEPATH"
    ml load htslib python/3.8.2
else
    ml load htslib python/3.6.1
fi
if [ "${avx2_flag}" -eq 2 ] ; then
    ml load plink2/20200314
else
    ml load plink2/20200314-non-AVX2
fi

pheno="/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.20200522.phe"
pheno_covar_col_end=94

if   [ "${pop}" == "related" ] ; then
    out="/oak/stanford/groups/mrivas/projects/related/"
    # out="/oak/stanford/groups/mrivas/projects/dev-related/"
    keep="/oak/stanford/groups/mrivas/ukbb24983/sqc/relatedness_20200514/semi_related.fam"
elif [ "${pop}" == "others" ] ; then
    out="/oak/stanford/groups/mrivas/projects/gwas_others/"
    # out="/oak/stanford/groups/mrivas/projects/dev-gwas_others/"
    keep="/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200522/ukb24983_others.phe"
else
    echo "unsupported population ${pop}" >&2 ; exit 1
fi

if [ "${regression_type}" == "linear" ] ; then
    start_idx=$( perl -e "print(int(((${batch_idx} - 1) * ${batch_size_qt}) + 1))" )
    end_idx=$(   perl -e "print(int(  ${batch_idx}      * ${batch_size_qt}))" )
    pheno_name=$(head -n1 ${pheno} | tr "\t" "\n" | egrep -n 'INI|QT_FC' | awk -v pheno_covar_col_end=${pheno_covar_col_end} -v FS=':' '(NR>pheno_covar_col_end){print $1}' | awk -v start_idx=${start_idx} -v end_idx=${end_idx} 'start_idx <= NR && NR <= end_idx' | tr "\n" "," | rev | cut -c2- | rev)
    plink_suffix="glm.linear"
elif [ "${regression_type}" == "logistic" ] ; then
    start_idx=$( perl -e "print(int(((${batch_idx} - 1) * ${batch_size_bin}) + 1))" )
    end_idx=$(   perl -e "print(int(  ${batch_idx}      * ${batch_size_bin}))" )
    pheno_name=$(head -n1 ${pheno} | tr "\t" "\n" | egrep -n -v 'INI|QT_FC' | awk -v pheno_covar_col_end=${pheno_covar_col_end} -v FS=':' '(NR>pheno_covar_col_end){print $1}' | awk -v start_idx=${start_idx} -v end_idx=${end_idx} 'start_idx <= NR && NR <= end_idx' | tr "\n" "," | rev | cut -c2- | rev)
    plink_suffix="glm.logistic.hybrid"
else
    echo "unsupported regression_type ${regression_type}" >&2 ; exit 1
fi

echo ${pheno_name} | tr "," '\n' | while read idx ; do
    GBE_ID=$(head -n1 ${pheno} | tr "\t" "\n" | awk -v idx=${idx} 'NR==idx')    
    out1=${out}/ukb24983_v2_hg19.${pheno_name}.array-combined.one_array.${GBE_ID}.${plink_suffix}
    out2=${out}/ukb24983_v2_hg19.${pheno_name}.array-combined.both_arrays.${GBE_ID}.${plink_suffix}
    out_gz=${out}/ukb24983_v2_hg19.${GBE_ID}.array-combined.${plink_suffix}.gz

    if [ -f ${out1} ] && [ -f ${out2} ] ; then
        {
            ! head -n1 ${out1}
            cat ${out1} ${out2} | grep -v '#' | sort --parallel 6 -k1,1V -k2,2n -k3,3
        } | bgzip -l9 -@6 > ${out_gz}
        if [ -f ${out_gz} ] ; then rm ${out1} ${out2} ; fi
        echo ${out_gz}
    fi
done

exit 0
# instructions
- There are 1569 QTs and 1899 binary phenotypes
- 156 array jobs for QTs, 475 array jobs for bins
- Yosuke --> all of QTs and binary 1-200
- Manny --> binary 201-475

bash 20200606-gwas-additional-pops-rename.sh related 1 linear
bash 20200606-gwas-additional-pops-rename.sh related 1
bash 20200606-gwas-additional-pops-rename.sh related 2
bash 20200606-gwas-additional-pops-rename.sh others 1 linear
bash 20200606-gwas-additional-pops-rename.sh others 1
bash 20200606-gwas-additional-pops-rename.sh others 2

{ for pop in related others ; do for idx in $(seq 2 156) ; do echo "[running] $pop linear $idx" ; bash 20200606-gwas-additional-pops-rename.sh $pop $idx linear ; done ; done ; for pop in related others ; do for idx in $(seq 3 200) ; do echo "[running] $pop logistic $idx" ; bash 20200606-gwas-additional-pops-rename.sh $pop $idx ; done ; done ; } | tee -a 20200606-rename.log

exit 0
# this rename was too slow. I started another set of commands in tmux sessions
for idx in $(seq 3 200) ; do echo "[running] others logistic $idx" ; bash 20200606-gwas-additional-pops-rename.sh others $idx ; done | tee -a 20200606-rename-others-logistic.log

for idx in $(seq 3 200) ; do echo "[running] related logistic $idx" ; bash 20200606-gwas-additional-pops-rename.sh related $idx ; done | tee -a 20200606-rename-related-logistic.log

for idx in $(seq 50 156) ; do echo "[running] others linear $idx" ; bash 20200606-gwas-additional-pops-rename.sh others $idx linear ; done | tee -a 20200606-rename-others-linear.log

###

[ytanigaw@sh02-02n07 /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas]$ find /oak/stanford/groups/mrivas/projects/related/ -name "ukb24983_v2_hg19.*-*.array-combined.*hybrid" | while read f ; do basename $f | awk -v FS='.' '{print $2}' ; done | sort6 -u
95-595
596-1096
1097-1597
1598-2098
2099-2599
2600-3100
# those files are failed sumstats --> delted



#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --output=logs/merge.%A.out
#SBATCH  --error=logs/merge.%A.err
#SBATCH --nodes=1
#SBATCH --cores=16
#SBATCH --mem=3000000
#SBATCH --time=0-20:00:00
#SBATCH -p mrivas,normal
set -beEuo pipefail

cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )

##############
tmp_dir_root=${LOCAL_SCRATCH}
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT
##############

cat_or_zstdcat () {
    local file=$1
    if [ "${file%.zst}.zst" == "${file}" ] ; then zstdcat $file ; else cat $file ; fi
}

prep_files () {
    local zstfile=$1
    local prefix=$2
    local tmp_dir=$3
    
    local pfile=$(echo $zstfile | sed -e "s/.zst$//g" | sed -e "s/.pvar$//g")
    local tmp_pvar="${tmp_dir}/${prefix}.tmp.in.pvar"
    local tmp_bim="${tmp_dir}/$(basename $pfile).bim"
    local tmp_bed="${tmp_dir}/$(basename $pfile).bed"
    local tmp_fam="${tmp_dir}/$(basename $pfile).fam"

    ln -s ${pfile}.fam ${tmp_fam}    
    ln -s ${pfile}.bed ${tmp_bed}
    
    cat_or_zstdcat ${zstfile} | grep -v '##' > ${tmp_pvar}
    Rscript assign_short_names_for_merge.R ${tmp_pvar} ${prefix} ${pfile}.IDs.tsv
    rm ${tmp_pvar}
    
    cat ${pfile}.IDs.tsv \
    | awk -v FS='\t' -v OFS='\t' '(NR>1){print $1, $6, 0, $2, $5, $4}' > ${tmp_bim}
}
##############

merge_list_tsv=$1
out="/oak/stanford/groups/mrivas/ukbb/24983/array_imp_combined/pgen/ukb24983_ukb24983_cal_hla_cnv_imp"

plink_opts="--memory ${mem} --threads ${cores}"
tmp_merge_list=${tmp_dir}/merge.lst

# prepare temp files for merge (assign short IDs)
cat ${merge_list_tsv} | while read f p ; do 
    echo "processing $f" >&2
    prep_files $f $p $tmp_dir
    echo ${tmp_dir}/$(basename $f) | sed -e "s/.zst$//g" | sed -e "s/.pvar$//g" >> ${tmp_merge_list}
done

cat ${tmp_merge_list} >&2

plink ${plink_opts} --merge-list ${tmp_merge_list} --make-bed --keep-allele-order --out ${out}

mv "${out}.bim" "${out}.shortnames.bim"

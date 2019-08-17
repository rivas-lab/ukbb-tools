#!/bin/bash
set -beEuo pipefail

# Given a GWAS summary statistics in PLINK format, 
# check the "REF" col against the reference sequence in a FASTA file 

__ref_hg19_fa='/oak/stanford/groups/mrivas/public_data/genomes/hg19/hg19.fa'

in_sumstats=$1
if [ $# -gt 1 ] ; then ref_fa=$2 ; else ref_fa="${__ref_hg19_fa}" ; fi

tmp_dir_root=/dev/shm/u/$USER
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT

cat_or_zcat () {
    local file=$1
    if [ "${file%.gz}.gz" == "${file}" ] ; then zcat ${file} ; else cat ${file} ; fi
}

show_header () {
    local file=$1
    ! cat_or_zcat $file | head -n1
}

check_flip_header () {
    local in_sumstats=$1
    echo "$(show_header $in_sumstats) FASTA_REF" | tr " " "\t"
}

get_col_idx () {
    local file=$1
    local key=$2
    show_header $file | tr "\t" "\n" | awk -v key=$key '($0 == key){print NR}'
}

check_flip_body () {
    local in_sumstats=$1
    local ref_fa=$2    
    local chr_prefix="chr"
    local field_sep="!"
    zcat ${in_sumstats} | tr "\t" "${field_sep}" \
        | awk -v OFS='\t' -v FS=${field_sep} -v chr="${chr_prefix}" \
        '(NR>1){print chr $1, $2-1, $2-1+length($4), $3, chr $0}' \
        | bedtools getfasta -fi ${ref_fa} -bed /dev/stdin -bedOut \
        | tr ${field_sep} "\t" | cut -f5- | sed -e "s/^${chr_prefix}//g"
}

check_flip () {
    local in_sumstats=$1
    local ref_fa=$2  
    local tmp_dir=$3

    local tmp_ref_fa=${tmp_dir}/$(basename $tmp_dir) 
    cp ${ref_fa}     ${tmp_ref_fa}
    cp ${ref_fa}.fai ${tmp_ref_fa}.fai
    date >&2
    check_flip_header ${in_sumstats}
    check_flip_body   ${in_sumstats} ${tmp_ref_fa}
    date >&2
}

fix_flip_binary () {
    local in_sumstats=$1
    local tmp_ref_fa=$2

    local col_REF=$(get_col_idx $in_sumstats "REF")
    local col_ALT=$(get_col_idx $in_sumstats "ALT")
    local col_A1=$(get_col_idx $in_sumstats "A1")
    local col_OR=$(get_col_idx $in_sumstats "OR")

    show_header ${in_sumstats}
    check_flip_body ${in_sumstats} ${tmp_ref_fa} \
        | awk -v FS='\t' -v OFS='\t' -v cOR=${col_OR} \
        -v cREF=${col_REF} -v cALT=${col_ALT} -v cA1=${col_A1} '
        (toupper($NF) == $cREF){ print $0 } 
        (toupper($NF) == $cA1 && toupper($NF) == $cALT){
            $cALT = $cREF ;
            $cREF = $cA1 ;
            $cA1  = $cALT ; 
            $cOR  = exp(-1 * log($cOR)) ;
            print $0
        }' \
        | rev | cut -f2- | rev
}

fix_flip_quantitative () {
    local in_sumstats=$1
    local tmp_ref_fa=$2

    local col_REF=$(get_col_idx $in_sumstats "REF")
    local col_ALT=$(get_col_idx $in_sumstats "ALT")
    local col_A1=$(get_col_idx $in_sumstats "A1")
    local col_BETA=$(get_col_idx $in_sumstats "BETA")

    show_header ${in_sumstats}
    check_flip_body ${in_sumstats} ${tmp_ref_fa} \
        | awk -v FS='\t' -v OFS='\t' -v cBETA=${col_BETA} \
        -v cREF=${col_REF} -v cALT=${col_ALT} -v cA1=${col_A1} '
        (toupper($NF) == $cREF){ print $0 } 
        (toupper($NF) == $cA1 && toupper($NF) == $cALT){
            $cALT  = $cREF ;
            $cREF  = $cA1 ;
            $cA1   = $cALT ; 
            $cBETA = -1 * $cBETA ;
            print $0
        }' \
        | rev | cut -f2- | rev
}

fix_flip () {
    local in_sumstats=$1
    local ref_fa=$2  
    local tmp_dir=$3

    local tmp_ref_fa=${tmp_dir}/$(basename $tmp_dir) 
    cp ${ref_fa}     ${tmp_ref_fa}
    cp ${ref_fa}.fai ${tmp_ref_fa}.fai

    local col_OR=$(  get_col_idx $in_sumstats "OR")
    local col_BETA=$(get_col_idx $in_sumstats "BETA")

    if   [ ! -z "${col_OR}" ] && [   -z "${col_BETA}" ] ; then
        fix_flip_binary ${in_sumstats} ${tmp_ref_fa}
    elif [   -z "${col_OR}" ] && [ ! -z "${col_BETA}" ] ; then
        fix_flip_quantitative ${in_sumstats} ${tmp_ref_fa}
    else
        echo "failed to detect the file type. Is this binary or quantitative?" &>2 ; 
        exit 1
    fi
}

#check_flip $in_sumstats $ref_fa ${tmp_dir}
fix_flip ${in_sumstats} ${ref_fa} ${tmp_dir}


#!/bin/bash
set -beEuo pipefail

[[ ${REF_FA_DIR:-}  -eq 1 ]]  && return || readonly REF_FA_DIR="/oak/stanford/groups/mrivas/public_data/genomes"
[[ ${REF_HG19_FA:-} -eq 1 ]]  && return || readonly REF_HG19_FA="${REF_FA_DIR}/hg19/hg19.fa"
[[ ${REF_HG38_FA:-} -eq 1 ]]  && return || readonly REF_HG38_FA="${REF_FA_DIR}/hg38/hg38.fa"

cat_or_zcat () {
    local file=$1
    if [ "${file%.gz}.gz" == "${file}" ] ; then 
        zcat ${file} 
    elif [ "${file%.zst}.zst" == "${file}" ] ; then 
        zstdcat ${file}
    else
        cat ${file}
    fi
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
    show_header $file | sed -e "s/^#//g" | tr "\t" "\n" | awk -v key=$key '($0 == key){print NR}'
}

check_flip_body () {
    local in_sumstats=$1
    local ref_fa=$2    
    local chr_prefix="chr"
    local field_sep="!"

    local col_CHROM=$(get_col_idx $in_sumstats "CHROM")
    local col_POS=$(get_col_idx $in_sumstats "POS")
    local col_ID=$(get_col_idx $in_sumstats "ID")
    local col_REF=$(get_col_idx $in_sumstats "REF")

    cat_or_zcat ${in_sumstats} | tr "\t" "${field_sep}" \
        | awk -v OFS='\t' -v FS=${field_sep} -v chr="${chr_prefix}" \
        -v cCHROM=${col_CHROM} -v cPOS=${col_POS} -v cID=${col_ID} -v cREF=${col_REF} \
        '(NR>1){print chr $cCHROM, $cPOS-1, $cPOS-1+length($cREF), $cID, chr $0}' \
        | sed -e "s/chrXY/chrX/g" \
        | sed -e "s/chrMT/chrM/g" \
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

    # we call check_flip_body to get the reference allele at the last column.
    # with that, we apply two awk scripts to apply flip-fix.
    # the first awk script ensures ALT == A1
    # the second awk script ensures REF is what you see in the fasta file

    show_header ${in_sumstats}
    check_flip_body ${in_sumstats} ${tmp_ref_fa} \
        | awk -v FS='\t' -v OFS='\t' \
        -v cREF=${col_REF} -v cALT=${col_ALT} -v cA1=${col_A1} '
        (toupper($cA1) == toupper($cALT)){ print $0 }
        (toupper($cA1) == toupper($cREF)){ 
            swap  = $cREF ;
            $cREF = $cALT ;
            $cALT = swap ;
            print $0 
        } ' \
        | awk -v FS='\t' -v OFS='\t' -v cOR=${col_OR} \
        -v cREF=${col_REF} -v cALT=${col_ALT} -v cA1=${col_A1} '
        (toupper($NF) == $cREF){ print $0 } 
        (toupper($NF) == $cA1 && toupper($NF) == $cALT){
            swap  = $cREF ;
            $cREF = $cA1 ;
            $cALT = swap ;
            $cA1  = swap ;
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
        | awk -v FS='\t' -v OFS='\t' \
        -v cREF=${col_REF} -v cALT=${col_ALT} -v cA1=${col_A1} '
        (toupper($cA1) == toupper($cALT)){ print $0 }
        (toupper($cA1) == toupper($cREF)){ 
            swap  = $cREF ;
            $cREF = $cALT ;
            $cALT = swap ;
            print $0 
        } ' \
        | awk -v FS='\t' -v OFS='\t' -v cBETA=${col_BETA} \
        -v cREF=${col_REF} -v cALT=${col_ALT} -v cA1=${col_A1} '
        (toupper($NF) == $cREF){ print $0 } 
        (toupper($NF) == $cA1 && toupper($NF) == $cALT){
            swap  = $cREF ;
            $cREF = $cA1 ;
            $cALT = swap ;
            $cA1  = swap ;
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
        echo "failed to detect the file type. Is this binary or quantitative?" >&2 ; 
        echo "col_OR=${col_OR} col_BETA=${col_BETA}"
        cat_or_zcat $in_sumstats | awk 'NR<5'
        exit 1
    fi
}


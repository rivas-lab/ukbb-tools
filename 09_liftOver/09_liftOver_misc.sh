#!/bin/bash
set -beEuo pipefail

# constants
[[ ${TO_BED_FIELD_SEP:-}  -eq 1 ]]  && return || readonly TO_BED_FIELD_SEP='!'
[[ ${BED_CHR_PREFIX:-}    -eq 1 ]]  && return || readonly BED_CHR_PREFIX='chr'

# reference data
# [[ ${REF_FA_DIR:-}  -eq 1 ]]  && return || readonly REF_FA_DIR="/oak/stanford/groups/mrivas/public_data/genomes"
[[ ${REF_FA_DIR:-}  -eq 1 ]]  && return || readonly REF_FA_DIR="/scratch/groups/mrivas/public_data/genomes"
[[ ${REF_HG19_FA:-} -eq 1 ]]  && return || readonly REF_HG19_FA="${REF_FA_DIR}/hg19/hg19.fa"
[[ ${REF_HG38_FA:-} -eq 1 ]]  && return || readonly REF_HG38_FA="${REF_FA_DIR}/hg38/hg38.fa"


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

get_ref_fa () {
    local assembly=$1
    if [ "${assembly}" == "hg19" ] || [ "${assembly}" == "grch37" ] ; then
        local ref_fa="${REF_HG19_FA}"
    elif [ "${assembly}" == "hg38" ] || [ "${assembly}" == "grch38" ] ; then
        local ref_fa="${REF_HG38_FA}"
    fi
    echo "${ref_fa}"
}

show_header () {
    local file=$1
    # ! cat_or_zcat $file | head -n1
    cat_or_zcat $file | egrep '^#'
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

convertToBed () {
    local in_file=$1
    local to_bed_field_sep=$2
    local bed_chr_prefix=$3

    local col_CHROM=$( get_col_idx $in_file "CHROM")
    local col_POS=$(   get_col_idx $in_file "POS")
    local col_ID=$(    get_col_idx $in_file "ID")
    local col_REF=$(   get_col_idx $in_file "REF")
    
    cat_or_zcat ${in_file} \
    | egrep -v '^#' \
    | tr "\t" "${to_bed_field_sep}" \
    | awk -v OFS='\t' -v FS=${to_bed_field_sep} -v chr="${bed_chr_prefix}" \
        -v cCHROM=${col_CHROM} -v cPOS=${col_POS} -v cID=${col_ID} -v cREF=${col_REF} \
        '{print chr $cCHROM, $cPOS-1, $cPOS-1+length($cREF), $cID, chr $0}' \
    | sed -e "s/^${bed_chr_prefix}${bed_chr_prefix}/${bed_chr_prefix}/g" \
    | sed -e 's/^chrXY/chrX/g' \
    | sed -e 's/^chrMT/chrM/g' \
    | sed -e 's/^chr23/chrX/g' \
    | sed -e 's/^chr24/chrY/g' \
    | sed -e 's/^chr25/chrX/g' \
    | sed -e 's/^chr26/chrM/g'
}

## 
#  Flip-fix related functions
##

check_flip_body () {
    local in_file=$1
    local ref_fa=$2    
    local to_bed_field_sep=$3
    local bed_chr_prefix=$4
    
    convertToBed ${in_file} ${to_bed_field_sep} ${bed_chr_prefix} \
    | bedtools getfasta -fi ${ref_fa} -bed /dev/stdin -bedOut \
    | tr ${to_bed_field_sep} "\t" | cut -f5- | sed -e "s/^${bed_chr_prefix}//g"
}

check_flip () {
    local in_sumstats=$1
    local ref_fa=$2
    local to_bed_field_sep=$3
    local bed_chr_prefix=$4

    check_flip_header ${in_sumstats}
    check_flip_body   ${in_sumstats} ${ref_fa} ${to_bed_field_sep} ${bed_chr_prefix}
}

fix_flip_binary () {
    local in_sumstats=$1
    local ref_fa=$2
    local to_bed_field_sep=$3
    local bed_chr_prefix=$4

    local col_REF=$(get_col_idx $in_sumstats "REF")
    local col_ALT=$(get_col_idx $in_sumstats "ALT")
    local col_A1=$(get_col_idx $in_sumstats "A1")
    local col_OR=$(get_col_idx $in_sumstats "OR")

    # we call check_flip_body to get the reference allele at the last column.
    # with that, we apply two awk scripts to apply flip-fix.
    # the first awk script ensures ALT == A1
    # the second awk script ensures REF is what you see in the fasta file

    check_flip_body ${in_sumstats} ${ref_fa} ${to_bed_field_sep} ${bed_chr_prefix} \
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
    local ref_fa=$2
    local to_bed_field_sep=$3
    local bed_chr_prefix=$4

    local col_REF=$(get_col_idx $in_sumstats "REF")
    local col_ALT=$(get_col_idx $in_sumstats "ALT")
    local col_A1=$(get_col_idx $in_sumstats "A1")
    local col_BETA=$(get_col_idx $in_sumstats "BETA")
    
    check_flip_body ${in_sumstats} ${ref_fa} ${to_bed_field_sep} ${bed_chr_prefix} \
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
    local to_bed_field_sep=$3
    local bed_chr_prefix=$4

    local col_OR=$(  get_col_idx $in_sumstats "OR")
    local col_BETA=$(get_col_idx $in_sumstats "BETA")

    if   [ ! -z "${col_OR}" ] && [   -z "${col_BETA}" ] ; then
        show_header ${in_sumstats}
        fix_flip_binary ${in_sumstats} ${ref_fa} ${to_bed_field_sep} ${bed_chr_prefix}
    elif [   -z "${col_OR}" ] && [ ! -z "${col_BETA}" ] ; then
        show_header ${in_sumstats}
        fix_flip_quantitative ${in_sumstats} ${ref_fa} ${to_bed_field_sep} ${bed_chr_prefix}
    else
        echo "failed to detect the file type. Is this binary or quantitative?" >&2 ; 
        echo "col_OR=${col_OR} col_BETA=${col_BETA}"
        cat_or_zcat $in_sumstats | awk 'NR<5'
        exit 1
    fi
}

liftOver_body () {
    local in_file=$1
    local src_genome=$2
    local dst_genome=$3
    local out_unmapped=$4
    local threads=$5
    local to_bed_field_sep=$6
    local bed_chr_prefix=$7

    local chain_file="${CHAIN_DIR}/${src_genome}To${dst_genome^}.over.chain.gz"
    local col_CHROM=$( get_col_idx $in_file "CHROM")
    local col_POS=$(   get_col_idx $in_file "POS")

    convertToBed ${in_file} ${to_bed_field_sep} ${bed_chr_prefix} \
    | liftOver /dev/stdin ${chain_file} /dev/stdout ${out_unmapped} \
    | tr ${to_bed_field_sep} '\t' \
    | awk -v FS='\t' -v OFS='\t' -v cCHROM=${col_CHROM} -v cPOS=${col_POS} \
        '{ $(cCHROM + 4) = $1 ; $(cPOS + 4) = $2 + 1 ; print $0 }' \
    | cut -f5- \
    | sed -e "s/^${bed_chr_prefix}//g" \
    | sort -k${col_CHROM},${col_CHROM}V -k${col_POS},${col_POS}n --parallel=${threads}

}

liftOverWrapper () {
    local in_file=$1
    local src_genome=$2
    local dst_genome=$3
    local out_mapped=$4
    local out_unmapped=$5
    local threads=$6
    local to_bed_field_sep=$7
    local bed_chr_prefix=$8

    local col_OR=$(  get_col_idx ${in_file} "OR")
    local col_BETA=$(get_col_idx ${in_file} "BETA")

    if   [ ! -z "${col_OR}" ] || [ ! -z "${col_BETA}" ] ; then
        local ref_fa=$(get_ref_fa "${dst_genome}")
        local tmp_intermediate="${out_mapped%.gz}.liftOver.gz"

        { show_header ${in_file}
         liftOver_body ${in_file} ${src_genome} ${dst_genome} ${out_unmapped%.gz} ${threads} ${to_bed_field_sep} ${bed_chr_prefix}  
        } | bgzip -l9 --threads ${threads} > ${tmp_intermediate}

        fix_flip ${tmp_intermediate} ${ref_fa} ${to_bed_field_sep} ${bed_chr_prefix}
    else
        show_header ${in_file}
        liftOver_body ${in_file} ${src_genome} ${dst_genome} ${out_unmapped%.gz} ${threads} ${to_bed_field_sep} ${bed_chr_prefix} 
    fi | bgzip -l9 --threads ${threads} > ${out_mapped%.gz}.gz

    gzip -9 ${out_unmapped%.gz}

    if   [ ! -z "${col_OR}" ] || [ ! -z "${col_BETA}" ] ; then
        rm ${tmp_intermediate}
    fi
}

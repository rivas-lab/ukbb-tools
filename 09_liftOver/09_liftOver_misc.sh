#!/bin/bash
set -beEuo pipefail

[[ ${TO_BED_FIELD_SEP:-}  -eq 1 ]]  && return || readonly TO_BED_FIELD_SEP='!'
[[ ${BED_CHR_PREFIX:-}    -eq 1 ]]  && return || readonly BED_CHR_PREFIX='chr'

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

convertToBed () {
    local in_file=$1
    if [ $# -gt 1 ] ; then 
        local to_bed_field_sep=$2
    else
        local to_bed_field_sep=${TO_BED_FIELD_SEP}
    fi
    if [ $# -gt 2 ] ; then 
        local bed_chr_prefix=$3
    else
        local bed_chr_prefix=${BED_CHR_PREFIX}
    fi

    local col_CHROM=$( get_col_idx $in_file "CHROM")
    local col_POS=$(   get_col_idx $in_file "POS")
    local col_ID=$(    get_col_idx $in_file "ID")
    local col_REF=$(   get_col_idx $in_file "REF")
    
    cat_or_zcat ${in_file} \
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
#     | egrep -v '^#' \

liftOverWrapper () {
    local in_file=$1
    local src_genome=$2
    local dst_genome=$3
    local out_mapped=$4
    local out_unmapped=$5
    local threads=$6
    if [ $# -gt 6 ] ; then local to_bed_field_sep=$7 ; else local to_bed_field_sep=${TO_BED_FIELD_SEP} ; fi
    if [ $# -gt 7 ] ; then local bed_chr_prefix=$8 ; else local bed_chr_prefix=${BED_CHR_PREFIX} ; fi

    local chain_file="${CHAIN_DIR}/${src_genome}To${dst_genome^}.over.chain.gz"
    local col_CHROM=$( get_col_idx $in_file "CHROM")
    local col_POS=$(   get_col_idx $in_file "POS")
        
    {
        # header line
        cat_or_zcat ${in_file} | egrep '^#'

        # body
        convertToBed ${in_file} ${to_bed_field_sep} ${bed_chr_prefix} \
        | liftOver /dev/stdin ${chain_file} /dev/stdout ${out_unmapped} \
        | tr ${to_bed_field_sep} '\t' \
        | awk -v FS='\t' -v OFS='\t' -v cCHROM=${col_CHROM} -v cPOS=${col_POS} \
            '(NR>1){ $(cCHROM + 4) = $1 ; $(cPOS + 4) = $2 + 1 ; print $0 }' \
        | cut -f5- \
        | sed -e "s/^${bed_chr_prefix}//g" \
        | sort -k${col_CHROM},${col_CHROM}V -k${col_POS},${col_POS}n --parallel=${threads}
    
    } | bgzip -l9 --threads ${threads} > ${out_mapped%.gz}.gz
}

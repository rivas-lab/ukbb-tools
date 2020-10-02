#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.0.1"
NUM_POS_ARGS="1"

############################################################
# input
############################################################

pvar_file=$1

############################################################
# functions
############################################################

vcf_header () {
cat <<- EOF
	##fileformat=VCFv4.3
	##fileDate=20200418
	##source=PLINKv2.00
	##INFO=<ID=ORIGINAL_VAR_ID,Number=1,Type=String,Description="The variant ID in the original data release">
	##contig=<ID=X,length=154930725>
	##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	DECOY
EOF
}

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


############################################################
# main
############################################################

ml load zstd

vcf_header
cat_or_zcat ${pvar_file} \
| grep -v '#' \
| awk -v OFS='\t' '{print "chr" $1, $2, $3, $4, $5, ".", ".", $6, "GT", "1"}'


# $ zstdcat /oak/stanford/groups/mrivas/ukbb24983/imp/pgen/ukb24983_im
# p_chrX_v3.pvar.zst | head
# ##INFO=<ID=ORIGINAL_VAR_ID,Number=1,Type=String,Description="The variant ID in the original data release">
# #CHROM  POS     ID      REF     ALT     INFO
# X       2699528 X:2699528:T:C   T       C       ORIGINAL_VAR_ID=rs750499943


# $ head test_vcf.vcf
# ##fileformat=VCFv4.3
# ##fileDate=20200418
# ##source=PLINKv2.00
# ##INFO=<ID=ORIGINAL_VAR_ID,Number=1,Type=String,Description="The variant ID in the original data release">
# ##contig=<ID=X,length=154930725>
# ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
# #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  2502845_2502845
# X       2699528 X:2699528:T:C   T       C       .       .       ORIGINAL_VAR_ID=rs750499943     GT      0


# zstdcat /oak/stanford/groups/mrivas/ukbb24983/imp/pgen/ukb24983_imp_chrX_v3.pvar.zst | grep -v '#' \
# | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, ".", ".", $6, "GT", "0"}' | head

# > | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, ".", ".", $6, "GT", "0"}' | head
# X       2699528 X:2699528:T:C   T       C       .       .       ORIGINAL_VAR_ID=rs750499943     GT      0
# X       2699555 X:2699555:C:A   C       A       .       .       ORIGINAL_VAR_ID=rs311165        GT      0
# X       2699559 X:2699559:A:G   A       G       .       .       ORIGINAL_VAR_ID=rs186009975     GT      0
# X       2699625 X:2699625:A:G   A       G       .       .       ORIGINAL_VAR_ID=rs6655038       GT      0

#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)

ml load python/2.7

pvar=/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_cALL_v2_hg19.pvar
vep_out=$(dirname $(dirname ${pvar}))/annotation_20201002/dev/$(basename ${pvar%.zst} .pvar).vep101-loftee
vep_in_vcf=${vep_out}.input.vcf
assembly=GRCh37

if [ ! -s ${vep_out}.tsv ] ; then

    # convert pvar to vcf
    if [ ! -d $(dirname ${vep_in_vcf}) ] ; then mkdir -p $(dirname ${vep_in_vcf}) ; fi

    ! bash $(dirname ${SRCDIR})/helpers/pvar.to.vcf.sh ${pvar} | head -n1007 \
        | sed -e 's/chrMT/chrM/g' | sed -e 's/chrXY/chrX/g' > ${vep_in_vcf}

    # run VEP w/ loftee
    bash $(dirname ${SRCDIR})/helpers/vep_loftee_wrapper.sh ${assembly} ${vep_in_vcf} ${vep_out}

    # clean-up
    mv ${vep_out} ${vep_out}.vcf
    rm ${vep_in_vcf}

    # tabulize the VEP/loftee VCF
    python /oak/stanford/groups/mrivas/software/loftee/src/tableize_vcf.py \
        --vcf ${vep_out}.vcf --output ${vep_out}.tsv \
        --vep_info Allele,Consequence,IMPACT,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,ALLELE_NUM,DISTANCE,STRAND,FLAGS,VARIANT_CLASS,SYMBOL_SOURCE,HGNC_ID,CANONICAL,MANE,TSL,APPRIS,CCDS,ENSP,SWISSPROT,TREMBL,UNIPARC,GENE_PHENO,SIFT,PolyPhen,DOMAINS,miRNA,HGVS_OFFSET,AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,AA_AF,EA_AF,gnomAD_AF,gnomAD_AFR_AF,gnomAD_AMR_AF,gnomAD_ASJ_AF,gnomAD_EAS_AF,gnomAD_FIN_AF,gnomAD_NFE_AF,gnomAD_OTH_AF,gnomAD_SAS_AF,MAX_AF,MAX_AF_POPS,CLIN_SIG,SOMATIC,PHENO,PUBMED,VAR_SYNONYMS,MOTIF_NAME,MOTIF_POS,HIGH_INF_POS,MOTIF_SCORE_CHANGE,TRANSCRIPTION_FACTORS

# Note: in version 2017 of the pipeline, we extracted the following fields: Existing_variation,Gene,Consequence,HGVSp,LoF,LoF_filter,LoF_flags,LoF_info,SYMBOL
# https://github.com/rivas-lab/ukbb-tools/blob/416442322f761978efb02f0b8aa9fd55b9d8a054/17_annotation/annotate_bims.sh#L82
#

fi


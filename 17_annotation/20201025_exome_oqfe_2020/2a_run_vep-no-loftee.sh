#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)

ml load python/2.7

idx=$1
idx_pad=$(perl -e "print(sprintf('%03d', ${idx}))")

assembly=GRCh38
# vep_in_vcf=/scratch/groups/mrivas/ukbb24983/exome-oqfe2020-annotation/input/split.UKBexomeOQFE.pvar.body.${idx_pad}.vcf
vep_in_vcf=/scratch/groups/mrivas/ukbb24983/exome/annotation/20201025_exome_oqfe_2020/input/split.UKBexomeOQFE.pvar.body.${idx_pad}.vcf
vep_out=$(dirname $(dirname ${vep_in_vcf}))/output_vep/$(basename ${vep_in_vcf%} .vcf).vep101-loftee

if [ ! -d $(dirname ${vep_out}) ] ; then mkdir $(dirname ${vep_out}) ; fi

if [ ! -s ${vep_out}.tsv ] ; then

    # run VEP w/ loftee
    bash /oak/stanford/groups/mrivas/users/${USER}/repos/rivas-lab/ukbb-tools/17_annotation/helpers/vep_loftee_wrapper.sh \
        --skip_loftee ${assembly} ${vep_in_vcf} ${vep_out}

    # clean-up
    mv ${vep_out} ${vep_out}.vcf

    # tabulize the VEP/loftee VCF
    python /oak/stanford/groups/mrivas/software/loftee/src/tableize_vcf.py \
        --vcf ${vep_out}.vcf --output ${vep_out}.tsv \
        --vep_info Allele,Consequence,IMPACT,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,ALLELE_NUM,DISTANCE,STRAND,FLAGS,VARIANT_CLASS,SYMBOL_SOURCE,HGNC_ID,CANONICAL,MANE,TSL,APPRIS,CCDS,ENSP,SWISSPROT,TREMBL,UNIPARC,GENE_PHENO,SIFT,PolyPhen,DOMAINS,miRNA,HGVS_OFFSET,AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,AA_AF,EA_AF,gnomAD_AF,gnomAD_AFR_AF,gnomAD_AMR_AF,gnomAD_ASJ_AF,gnomAD_EAS_AF,gnomAD_FIN_AF,gnomAD_NFE_AF,gnomAD_OTH_AF,gnomAD_SAS_AF,MAX_AF,MAX_AF_POPS,CLIN_SIG,SOMATIC,PHENO,PUBMED,VAR_SYNONYMS,MOTIF_NAME,MOTIF_POS,HIGH_INF_POS,MOTIF_SCORE_CHANGE,TRANSCRIPTION_FACTORS
#        --vep_info Allele,Consequence,IMPACT,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,ALLELE_NUM,DISTANCE,STRAND,FLAGS,VARIANT_CLASS,SYMBOL_SOURCE,HGNC_ID,CANONICAL,MANE,TSL,APPRIS,CCDS,ENSP,SWISSPROT,TREMBL,UNIPARC,GENE_PHENO,SIFT,PolyPhen,DOMAINS,miRNA,HGVS_OFFSET,AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,AA_AF,EA_AF,gnomAD_AF,gnomAD_AFR_AF,gnomAD_AMR_AF,gnomAD_ASJ_AF,gnomAD_EAS_AF,gnomAD_FIN_AF,gnomAD_NFE_AF,gnomAD_OTH_AF,gnomAD_SAS_AF,MAX_AF,MAX_AF_POPS,CLIN_SIG,SOMATIC,PHENO,PUBMED,VAR_SYNONYMS,MOTIF_NAME,MOTIF_POS,HIGH_INF_POS,MOTIF_SCORE_CHANGE,TRANSCRIPTION_FACTORS,LoF,LoF_filter,LoF_flags,LoF_info

# Note: in version 2017 of the pipeline, we extracted the following fields: Existing_variation,Gene,Consequence,HGVSp,LoF,LoF_filter,LoF_flags,LoF_info,SYMBOL
# https://github.com/rivas-lab/ukbb-tools/blob/416442322f761978efb02f0b8aa9fd55b9d8a054/17_annotation/annotate_bims.sh#L82
#
fi

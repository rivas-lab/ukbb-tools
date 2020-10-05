#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)

ml load python/2.7 R/3.6 gcc

idx=$1
idx_pad=$(perl -e "print(sprintf('%03d', ${idx}))")

assembly=GRCh37

split_pvar=/scratch/groups/mrivas/ukbb24983/cnv/annotation_20201003/input/split.cnv.pvar.body.${idx_pad}.pvar
split_seq_pvar=${split_pvar%.pvar}.seq.pvar
vep_in_vcf=${split_pvar%.pvar}.vcf
vep_out=$(dirname $(dirname ${vep_in_vcf}))/output_vep/$(basename ${vep_in_vcf%} .vcf).vep101-loftee

# reference data
public_d="/oak/stanford/groups/mrivas/public_data"
vep_data="${public_d}/vep/20200912"
loftee_data=$(dirname ${vep_data})/20201002_loftee_data


if [ ! -s ${vep_out}.tsv ] ; then
    # add sequence to pvar
    bash cnv_pvar_fill_seq.sh ${split_pvar} ${split_seq_pvar}

    # convert pvar to vcf
    bash /oak/stanford/groups/mrivas/users/${USER}/repos/rivas-lab/ukbb-tools/17_annotation/helpers/pvar.to.vcf.sh ${split_seq_pvar} \
    | sed -e 's/chrMT/chrM/g' | sed -e 's/chrXY/chrX/g' > ${vep_in_vcf}

    # run VEP w/ loftee
    bash /oak/stanford/groups/mrivas/users/${USER}/repos/rivas-lab/ukbb-tools/17_annotation/helpers/vep_loftee_wrapper.sh \
    --vep_data ${vep_data} --loftee_data ${loftee_data} ${assembly} ${vep_in_vcf} ${vep_out}

    # clean-up
    mv ${vep_out} ${vep_out}.vcf

    # tabulize the VEP/loftee VCF
    python /oak/stanford/groups/mrivas/software/loftee/src/tableize_vcf.py \
        --vcf ${vep_out}.vcf --output ${vep_out}.tmp.tsv \
        --include_id --vep_info Allele,Consequence,IMPACT,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,ALLELE_NUM,DISTANCE,STRAND,FLAGS,VARIANT_CLASS,SYMBOL_SOURCE,HGNC_ID,CANONICAL,MANE,TSL,APPRIS,CCDS,ENSP,SWISSPROT,TREMBL,UNIPARC,GENE_PHENO,SIFT,PolyPhen,DOMAINS,miRNA,HGVS_OFFSET,AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,AA_AF,EA_AF,gnomAD_AF,gnomAD_AFR_AF,gnomAD_AMR_AF,gnomAD_ASJ_AF,gnomAD_EAS_AF,gnomAD_FIN_AF,gnomAD_NFE_AF,gnomAD_OTH_AF,gnomAD_SAS_AF,MAX_AF,MAX_AF_POPS,CLIN_SIG,SOMATIC,PHENO,PUBMED,VAR_SYNONYMS,MOTIF_NAME,MOTIF_POS,HIGH_INF_POS,MOTIF_SCORE_CHANGE,TRANSCRIPTION_FACTORS,LoF,LoF_filter,LoF_flags,LoF_info

# Note: in version 2017 of the pipeline, we extracted the following fields: Existing_variation,Gene,Consequence,HGVSp,LoF,LoF_filter,LoF_flags,LoF_info,SYMBOL
# https://github.com/rivas-lab/ukbb-tools/blob/416442322f761978efb02f0b8aa9fd55b9d8a054/17_annotation/annotate_bims.sh#L82
#
    
    mv ${vep_out}.tmp.tsv.log ${vep_out}.tsv.log

    # Allele strings are too long for CNV dataset. Let's clean them up.
    Rscript /dev/stdin ${vep_out}.tmp.tsv ${vep_out}.tsv << EOF
    suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))
    args <- commandArgs(trailingOnly=TRUE)
    in_f <- args[1]
    out_f <- args[2]

    in_f %>% fread() %>% select(-REF, -ALT) %>%
    mutate(Allele = if_else(str_length(Allele)>1, '+', '-')) %>%
    mutate(HGVSc = str_replace_all(HGVSc, 'ins[ACGT]+', 'ins')) %>%
    fwrite(out_f, sep='\t', na = "NA", quote=F)
EOF

    rm ${vep_out}.tmp.tsv
    bgzip ${vep_out}.vcf
fi

#!/bin/bash

# Error message
usage () {
    echo "$0: Annotates BIM-style files with VEP annotations in the GBE format."
    echo "Note: The input files must have at least 5 columns - one each for CHROM,"
    echo "POS, REF, ALT, and MAF, and also have a single-dot extension"
    echo "(e.g. .txt, .tsv, and not .txt.tsv)."
    echo "Below is the top of an example file:: "
    echo ""
    echo "CHROM   POS     REF     ALT    MAF"
    echo "1      100000051       A       G    3.1000e-2"
    echo ""
    echo "Minimal usage: sbatch $0 <BIM input filename> <assembly name (GRCh37 or GRCh38)>"
    echo "Suggested usage: sbatch --mem=64000 -t 1-00:00:00 -p mrivas --nodes=1 --cores=8 $0 <BIM input filename> <assembly name (GRCh37 or GRCh38)>"
    echo "For reference, annotating ~15M variants with these parameters takes ~7 hours."
}

if [ $# -ne 2 ] ; then usage>&2 ; exit 1 ; fi

echo "Loading requisite packages..."
ml load python/2.7
ml load biology
ml load vcftools
ml load samtools
ml load bcftools
ml load perl

input=$1
echo "Input:" $input

assembly_in=$2
echo "Assembly:" $assembly_in

temp_vcf="${input%.*}"_temp.vcf
vcf="${input%.*}".vcf.gz

output="${input%.*}"_vep.vcf
echo "Output:" $output
stats="${input%.*}"_vep.html

echo "Converting to VCF..."
python convert_to_vcf.py $input $temp_vcf
cat $temp_vcf | vcf-sort | bgzip -c > $vcf

echo "Removing temp file..."
rm $temp_vcf

echo "Tabixing..."
tabix -p vcf $vcf

echo "Running VEP..."
perl /oak/stanford/groups/mrivas/users/guhan/software/ensembl-vep/vep.pl \
        --input_file $vcf \
        --force \
        --fork 4 \
        --vcf \
        --cache \
        --dir_cache /oak/stanford/groups/mrivas/public_data/vep_cache_20170410 \
        --offline \
        --database \
        --allele_number \
        --fasta /oak/stanford/groups/mrivas/public_data/vep_cache_20170410/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
        --everything \
        --assembly $assembly_in \
        --plugin LoF,loftee_path:/oak/stanford/groups/mrivas/public_data/loftee,human_ancestor_fa:/oak/stanford/groups/mrivas/public_data/loftee_human_ancestor_20170411/human_ancestor.fa.gz,conservation_file:/oak/stanford/groups/mrivas/public_data/loftee/phylocsf_gerp.sql \
        --dir_plugins /oak/stanford/groups/mrivas/public_data/loftee \
        --output_file $output \
        --stats_file $stats
echo "Result written to" $output"."

#rm $vcf*

echo "Zipping and tabixing output..."
bgzip -f $output
tabix -p vcf $output.gz

echo "Exporting as .tsv..."
python /oak/stanford/groups/mrivas/public_data/loftee/src/tableize_vcf.py \
        --vcf $output.gz \
        --output "${input%.*}"_vep.tsv \
        --vep_info Existing_variation,Gene,Consequence,HGVSp,LoF,LoF_filter,LoF_flags,LoF_info,SYMBOL

zcat $output.gz | tail -n +2 | gzip > $output.gz.tmp && mv $output.gz.tmp $output.gz

echo "Adding consequence and maf fields to .tsv..."
python add_fields.py "${input%.*}"_vep.tsv $output.gz $input "${input%.*}"_cf_vep.tsv

echo "Shuffling and adding columns..."
# Renoves filter field, shuffles gene symbol to later, then maf
awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13,$15,$14,$16}' "${input%.*}"_cf_vep.tsv > "${input%.*}"_cf_vep.tsv_tmp

# Rename some columns
sed -i '1s/Existing_variation/ID/g' "${input%.*}"_cf_vep.tsv_tmp && mv "${input%.*}"_cf_vep.tsv_tmp "${input%.*}"_cf_vep.tsv
sed -i '1s/SYMBOL/Gene_symbol/g' "${input%.*}"_cf_vep.tsv
sed -i '1s/MAF/maf/g' "${input%.*}"_cf_vep.tsv

# Join on gnomad AF


# Add additional fields - figure out gnomad_af
#awk -F'\t' 'BEGIN{OFS="\t"}{if (NR==1) {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,"f_miss","f_miss_bileve","f_miss_wcsg","freq","hwe_p",$15,"ld_indep","wcsg_only","bileve_only","filter","missingness","hwe","mcpi","gnomad_af","mgi","mgi_notes","all_filters","Gene_symbol"} else {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,"NA","NA","NA",$15,"NA",$15,"NA","NA","NA","0","0","0","NA","NA","NA","NA","0",$14}}' "${input%.*}"_cf_vep.tsv > "${input%.*}"_variant_annots.tsv

# Replace NA with quotes
#sed -i 's/\<NA\>/""/g' "${input%.*}"_variant_annots.tsv

#echo "Compressing..."
#gzip -f "${input%.*}"_variant_annots.tsv

#rm "${input%.*}"_cf_vep.tsv "${input%.*}"_vep*

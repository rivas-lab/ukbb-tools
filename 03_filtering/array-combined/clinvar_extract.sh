#!/bin/bash
set -beEuo pipefail

ml load bcftools

clinvar_vcf='/oak/stanford/groups/mrivas/public_data/ClinVar/vcf_GRCh37/clinvar_20200914.vcf.gz'

# this script extracts the Pathogenic and likely pathogenic variants from ClinVar

{

echo "#CHROM POS REF ALT CLNSIG"

bcftools query \
  -i '(FILTER == ".") && (INFO/CLNSIG ~ "pathogenic/i") && (INFO/CLNSIG !~ "Conflicting_interpretations_of_pathogenicity/i")' \
  -f '%CHROM %POS %REF %ALT %CLNSIG\n' ${clinvar_vcf} ${clinvar_vcf%.vcf.gz}_papu.vcf.gz \
  | grep -v 'Conflicting_interpretations_of_pathogenicity'

} | tr ' ' '\t' > clinvar_20200914_patho.tsv

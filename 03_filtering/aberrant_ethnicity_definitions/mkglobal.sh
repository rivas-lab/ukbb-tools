## Generate each self identity group

##################################################################
# 1) Derive the self-reported ethnic group phenotype definitions #
##################################################################
# If you are applying to a different application, then you can extract the relevant column directly
#cut -f 1,3986 /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/9796/24611/download/ukb24611.tab > ethnicity.phe

# In Rivas lab, there's already a derived file with the phenotypes
cut -f 1-2 /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/self_reported_ethnicity/misc/ukb9796_ukb24611_f21000.tsv > ethnicity.phe

# List of all individual IDs in the application (for creating exclusion list of "everyone else")
awk '{print $1 "\t" $1;}' ethnicity.phe > all_individuals.phe

WHOLE_GENOME_GENOTYPES_PREFIX=/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_cALL_v2

##################################################################
# 2) extract individuals of each self-identity group of interest #
##################################################################
# Remove White British individuals because they have already been defined
grep '	1'  ethnicity.phe | grep -v '1001$' | cut -f 1-2 > self_identified_white.txt
grep '	1'  ethnicity.phe | grep -v '1001$' | cut -f 1 | awk '{print $1 "\t" $1;}' > self_identified_white.phe

grep '	3'  ethnicity.phe | cut -f 1-2 > self_identified_southasian.txt
grep '	3'  ethnicity.phe | cut -f 1 | awk '{print $1 "\t" $1;}' > self_identified_southasian.phe

grep '	4'  ethnicity.phe | cut -f 1-2 > self_identified_african.txt
grep '	4'  ethnicity.phe | cut -f 1 | awk '{print $1 "\t" $1;}' > self_identified_african.phe

## East Asian individuals are in two places in the self identity
grep '	5'  ethnicity.phe | cut -f 1-2 > self_identified_chinese.txt
grep '	5'  ethnicity.phe | cut -f 1 | cut -f 1 | awk '{print $1 "\t" $1;}' > self_identified_chinese.phe

# Optional, if you want to include the non-Chinese-identified individuals with ancestry that is similar to that of the Chinese
# and who self-identify as "Any other Asian background" -- we do not apply this.
#grep '	3004'  ethnicity.phe | cut -f 1-2 >> self_identified_chinese.txt
#grep '	3004'  ethnicity.phe | cut -f 1 | cut -f 1 | awk '{print $1 "\t" $1;}' >> self_identified_chinese.phe

grep '	200[12]'  ethnicity.phe | cut -f 1-2 > self_identified_mixedwhiteafrican.txt
grep '	200[12]'  ethnicity.phe | cut -f 1 | cut -f 1 | awk '{print $1 "\t" $1;}' > self_identified_mixedwhiteafrican.phe

grep '	2003'  ethnicity.phe | cut -f 1-2 > self_identified_mixedwhiteasian.txt
grep '	2003'  ethnicity.phe | cut -f 1 | cut -f 1 | awk '{print $1 "\t" $1;}' > self_identified_mixedwhiteasian.phe

######################################################################
# 3) pre-process the individuals to be QC'd out, as well as PCA-able #
######################################################################
module load R/3.4.0
# not_used_in_pca_calculation.phe list
Rscript filter_indiv_pca.R
# pass_qc.phe list
Rscript filter_indiv_exclude.R

# load a very compatible version of plink2
module load plink2/20190402-non-AVX2

#####################################################################################
# 4) run a pre-PCA of self-identified individuals to define major axes of variation #
#####################################################################################

## Run pre-exact PCA for small sample sizes
ls self_identified_{southasian,african,chinese,mixedwhiteafrican,mixedwhiteasian}.phe | while read selfidentity; do
    plink2 --pca var-wts 10 vcols=chrom,pos,ref,alt1,alt,maj,nonmaj --bfile ${WHOLE_GENOME_GENOTYPES_PREFIX} --keep $selfidentity --out `basename $selfidentity .phe`.firstpca10 --maf 0.1 --bp-space 10000 --remove not_used_in_pca_calculation.phe --seed 42 --threads 16
done

## Run pre-approx PCA for large sample sizes
ls self_identified_white.phe | while read selfidentity; do
    plink2 --pca approx var-wts 10 vcols=chrom,pos,ref,alt1,alt,maj,nonmaj --bfile ${WHOLE_GENOME_GENOTYPES_PREFIX} --keep $selfidentity --out `basename $selfidentity .phe`.firstpca10 --maf 0.1 --bp-space 10000 --remove not_used_in_pca_calculation.phe --seed 42 --threads 16
done

# project all samples into the PCs for each ethnic group
ls self_identified_{white,southasian,african,chinese,mixedwhiteafrican,mixedwhiteasian}.phe | while read selfidentity; do
    plink2 --bfile ${WHOLE_GENOME_GENOTYPES_PREFIX} --score `basename $selfidentity .phe`.firstpca10.eigenvec.var 3 6 header-read no-mean-imputation --score-col-nums 9-18 --out `basename $selfidentity .phe`.firstpca10.projection --threads 16
    cat `basename $selfidentity .phe`.firstpca10.projection.sscore | tr -d '#' | cut -f 1-2,5- > `basename $selfidentity .phe`.firstpca10.phe
done

###########################################################################################
# 5) run aberrant to detect outliers in the set of self-identified individuals in pre-PCA #
###########################################################################################
mkdir lambda10
## Filter using the PC1/PC2, PC3/PC4, PC5/PC6 strategy used in the UKBB marker paper
Rscript aberrant_self_identity.R

####################################################################################################
# 6) run post-PCA to define axes of variation relevant to non-outlier, self-identified individuals #
####################################################################################################
# Run exact PCA for small sample sizes
ls lambda10/inlier.self_identified_{african,chinese,mixedwhiteafrican,southasian,mixedwhiteasian}.phe | while read phe; do
    plink2 --pca var-wts 20 vcols=chrom,pos,ref,alt1,alt,maj,nonmaj --bfile ${WHOLE_GENOME_GENOTYPES_PREFIX} --keep $phe --out `dirname $phe`/`basename $phe .phe`.pca20 --maf 0.1 --bp-space 10000 --remove not_used_in_pca_calculation.phe --seed 42 --threads 16
done

# Run approx PCA for large sample size self-identified White individuals
ls lambda10/inlier.self_identified_white.phe | while read phe; do
    plink2 --pca approx var-wts 20 vcols=chrom,pos,ref,alt1,alt,maj,nonmaj --bfile ${WHOLE_GENOME_GENOTYPES_PREFIX} --keep $phe --out `dirname $phe`/`basename $phe .phe`.pca20 --maf 0.1 --bp-space 10000 --remove not_used_in_pca_calculation.phe --seed 42 --threads 16
done

# project all samples into the PCs for each ethnicity
ls lambda10/inlier.self_identified_{white,african,chinese,mixedwhiteafrican,southasian,mixedwhiteasian}.phe | while read phe; do
    plink2 --bfile ${WHOLE_GENOME_GENOTYPES_PREFIX} --score `dirname $phe`/`basename $phe .phe`.pca20.eigenvec.var 3 6 header-read no-mean-imputation --score-col-nums 9-28 --out `dirname $phe`/`basename $phe .phe`.pca20.projection --threads 16
    cat `dirname $phe`/`basename $phe .phe`.pca20.projection.sscore | tr -d '#' | cut -f 1-2,5- > `dirname $phe`/`basename $phe .phe`.pca20.phe
done

#####################################################################################################
# 6) run KING to get a list of individuals who are unrelated in each of the aberrant-defined groups #
#####################################################################################################
ls lambda10/inlier.self_identified_{african,chinese,mixedwhiteafrican,southasian,mixedwhiteasian}.phe | while read indiv; do
    plink2 --king-cutoff 0.04419417 --bfile ${WHOLE_GENOME_GENOTYPES_PREFIX} --geno 0.2 --out `dirname $indiv`/`basename $indiv .phe`.king --threads 16 --maf 0.1 --hwe 1e-6 midp --seed 42 --keep $indiv
done

# filter to individuals who pass QC
ls lambda10/*.in.id | while read id; do
    grep -Ff pass_qc.phe $id > `dirname $id`/`basename $id .id`.qc.id
done

# flip file to produce one suitable for --exclude (in addition to the *.qc.id for --keep)
ls lambda10/*.qc.id | while read king; do
    grep -vFf $king all_individuals.phe > `dirname $king`/`basename $king .id`.exclude.id
done

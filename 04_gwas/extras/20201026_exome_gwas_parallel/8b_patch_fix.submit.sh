#!/bin/bash
set -beEuo pipefail

# wc_l_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l.20201104-090826.tsv
# wc_l_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l.20201104-121324.tsv
# wc_l_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l.20201104-133054.tsv
# wc_l_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l.20201104-171349.tsv
wc_l_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l.20201104-175636.tsv

fs=(
/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/white_british/ukb24983_exomeOQFE.HC1275.glm.logistic.hybrid.gz 
)

# cat ${wc_l_f} | awk -v FS='\t' '(NR>1&&$4 != 1){print $NF}' \
cat ${wc_l_f} | awk -v FS='\t' '( (NR>1) && (! ($4 == 1 && ($5 == 17777817 || $5 == 17777951)) ) ){print $NF}' \
| egrep -v $(for f in ${fs[@]} ; do echo $f ; done | tr '\n' '|' | rev | cut -c2- | rev) \
| while read f ; do

echo $f

sbatch -p mrivas,normal --nodes=1 --mem=80000 --cores=4 --time=2:00:00 \
--job-name=fix --output=logs/fix.%A.out --error=logs/fix.%A.err \
8_patch_fix.wrapper.sh \
${f} 4 80000

done

exit 0

###################################

# we submitted 24cpus, 64GB jobs for the following 6 files

fs=(
/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/others/ukb24983_exomeOQFE.HC1257.glm.logistic.hybrid.gz
/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/white_british/ukb24983_exomeOQFE.BIN_FC10003591.glm.logistic.hybrid.gz
/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/white_british/ukb24983_exomeOQFE.BIN_FC40001508.glm.logistic.hybrid.gz
/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/white_british/ukb24983_exomeOQFE.HC869.glm.logistic.hybrid.gz
/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/white_british/ukb24983_exomeOQFE.BIN22127.glm.logistic.hybrid.gz
/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/white_british/ukb24983_exomeOQFE.HC1158.glm.logistic.hybrid.gz
/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/white_british/ukb24983_exomeOQFE.HC1275.glm.logistic.hybrid.gz
)

for f in ${fs[@]} ; do
sbatch -p mrivas --qos=high_p --nodes=1 --mem=64000 --cores=24 --time=7-0:00:00 \
--job-name=fix --output=logs/fix.%A.out --error=logs/fix.%A.err \
8_patch_fix.wrapper.sh \
${f} 24 60000
done

fs=(
/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/white_british/ukb24983_exomeOQFE.HC1275.glm.logistic.hybrid.gz
)

for f in ${fs[@]} ; do
sbatch -p mrivas --qos=high_p --nodes=1 --mem=64000 --cores=24 --time=7-0:00:00 \
--job-name=fix --output=logs/fix.%A.out --error=logs/fix.%A.err \
8_patch_fix.wrapper.WB_HC1275.sh \
${f} 24 60000
done


fs=(
/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/related/ukb24983_exomeOQFE.BIN_FC4006150.glm.logistic.hybrid.gz
)

for f in ${fs[@]} ; do
sbatch -p mrivas,normal,owners --nodes=1 --mem=40000 --cores=4 --time=2:00:00 \
--job-name=fix --output=logs/fix.%A.out --error=logs/fix.%A.err \
8_patch_fix.wrapper.sh \
${f} 4 40000
done

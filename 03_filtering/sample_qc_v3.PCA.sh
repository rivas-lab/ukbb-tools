#!/bin/bash
set -beEuo pipefail

in_pfile="/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_cALL_v2_hg19"
pop_str_d="/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20190805"

mem=120000
cpu=10

pops=("e_asian" "s_asian" "african" "non_british_white")
out_d="${pop_str_d}/pca"
if [ ! -d ${out_d} ] ; then mkdir -p ${out_d} ; fi

# create a temp directory
tmp_dir_root=${LOCAL_SCRATCH}
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT

tmp_pvar=${tmp_dir}/var_filter
tmp_prune=${tmp_dir}/indep
tmp_bed=${tmp_dir}/pca_in
tmp_pca_d=${tmp_dir}/pca
mkdir -p ${tmp_pca_d}

plink_common_opts=" --memory ${mem} --threads ${cpu}"

# variant filter
echo "6 25477797 36448354 MHC" | tr " " "\t" \
    | plink2 ${plink_common_opts} --geno 0.1 --var-filter --snps-only just-acgt \
    --chr 1-22 --maf 0.05 --max-maf 0.95 --hwe 1e-10 midp --max-alleles 2 \
    --rm-dup exclude-all --exclude bed1 /dev/stdin \
    --make-just-pvar zs --pfile ${in_pfile} --out ${tmp_pvar}

# LD pruning
zstdgrep -v '#' ${tmp_pvar}.pvar.zst | cut -f3 \
    | plink2 ${plink_common_opts} --indep-pairwise 50 5 .5 --extract /dev/stdin \
    --out ${tmp_prune} --pfile ${in_pfile}

for pop in ${pops[@]} ; do
keep_f="${pop_str_d}/ukb24983_${pop}.phe"
plink2 ${plink_common_opts} --pfile ${in_pfile} --extract ${tmp_prune}.prune.in \
    --pca 10 var-wts approx vzs vcols=chrom,pos,ref,alt1,alt,maj,nonmaj \
    --keep ${keep_f} --out ${out_d}/ukb24983_${pop}_pca --seed 20190805
cat ${tmp_pvar}.log ${tmp_prune}.log ${out_d}/ukb24983_${pop}_pca.log  > ${out_d}/ukb24983_${pop}_pca.plink.log
rm ${out_d}/ukb24983_${pop}_pca.log
done

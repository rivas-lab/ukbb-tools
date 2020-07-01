#!/bin/bash
set -beEuo pipefail

gz_file_symlink=$1
out_d="gwas-current-gz-wc"

out="${out_d}/$(basename $(dirname $gz_file_symlink)).$(basename $gz_file_symlink .gz).txt"

if [ ! -f ${out} ] ; then
    gz_file=$(readlink -f $gz_file_symlink)
    wc_l=$(zcat $gz_file | wc -l)
    echo $gz_file_symlink $gz_file $wc_l | tr ' ' '\t' > ${out}
fi
exit 0
####################################

cd /oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current
find  $(pwd) -name "*.gz" -type l | sort > /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/extras/202006-GWAS-finish/gwas-current-gz-list.$(date +%Y%m%d-%H%M%S).txt

There are 22073 files.

sbatch -p mrivas,owners,normal --nodes=1 --mem=4000 --cores=1 --time=30:00 --job-name=wc --output=logs/wc.%A_%a.out --error=logs/wc.%A_%a.err --array=1-960 /oak/stanford/groups/mrivas/users/ytanigaw/repos/yk-tanigawa/resbatch/parallel-sbatch.sh gwas-current-gz-wc.sh gwas-current-gz-list.20200627-190224.txt 23
Submitted batch job 3092539

grep 'array-end' logs/wc*.err  | awk -v FS='=' '{print $NF}' > gwas-current-gz-wc.finished.$(date +%Y%m%d-%H%M%S).lst

cat gwas-current-gz-wc.finished.20200629-100602.lst  | sort | comm -3 <(seq 1 960 | sort) /dev/stdin | tr "\n" "," | rev | cut -c2- | rev

sbatch -p mrivas,normal --nodes=1 --mem=4000 --cores=1 --time=30:00 --job-name=wc --output=logs/wc.%A_%a.out --error=logs/wc.%A_%a.err --array=14,15,16,17,18,19,20,21,88,89,90,91,92,93,94,95 /oak/stanford/groups/mrivas/users/ytanigaw/repos/yk-tanigawa/resbatch/parallel-sbatch.sh gwas-current-gz-wc.sh gwas-current-gz-list.20200627-190224.txt 23


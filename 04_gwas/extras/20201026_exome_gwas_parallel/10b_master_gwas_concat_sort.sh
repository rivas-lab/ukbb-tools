#!/bin/bash
set -beEuo pipefail

pop=$1
if [ $# -gt 1 ] ; then chrom=$2 ; else chrom="ALL" ; fi

ml load R/3.6 gcc

data_d="/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE"

if [ "${chrom}" == "ALL" ] ; then
    out_f=${data_d}/ukb24983_exomeOQFE.${pop}.p1e-3.tsv.gz
else
    out_f=${data_d}/ukb24983_exomeOQFE.${pop}.chr${chrom}.p1e-3.tsv.gz
fi

if [ ! -f ${out_f} ] ; then

    if [ "${pop}" == "metal" ] ; then

        TMPDIR=$SCRATCH/tmp/Rtmp Rscript 10_master_gwas_concat_sort.R \
        ${out_f%.gz} \
        ${chrom} \
        $(dirname ${data_d})/$(basename ${data_d})_filtered_p1e-3/${pop}.tsv.gz \

    else

        TMPDIR=$SCRATCH/tmp/Rtmp Rscript 10_master_gwas_concat_sort.R \
        ${out_f%.gz} \
        ${chrom} \
        $(dirname ${data_d})/$(basename ${data_d})_filtered_p1e-3/${pop}.linear.tsv.gz \
        $(dirname ${data_d})/$(basename ${data_d})_filtered_p1e-3/${pop}.logistic.hybrid.tsv.gz

    fi
    bgzip -l9 -@6 ${out_f%.gz}
fi

if [ ! -f ${out_f}.tbi ] ; then
    tabix -s1 -b2 -e2 -c'#' ${out_f}
fi
echo ${out_f}

exit 0
#######################
bash 10b_master_gwas_concat_sort.sh african

for pop in 'white_british' 'non_british_white' 's_asian' 'e_asian' 'related' 'others' 'metal' ; do
sbatch -p mrivas --qos=high_p --nodes=1 --mem=64000 --cores=6 --time=7:00:00 \
    --job-name=concat_sort.${pop} --output=logs/concat_sort.%A.out --error=logs/concat_sort.%A.err \
    10b_master_gwas_concat_sort.sh ${pop}
done

WB and metal --> out of memory error.

for pop in 'white_british' 'metal' ; do
for chr in Y X $(seq 1 22) ; do
sbatch -p mrivas --qos=high_p --nodes=1 --mem=64000 --cores=6 --time=7:00:00 \
    --job-name=concat_sort.${pop}.${chr} --output=logs/concat_sort.${pop}.${chr}.%A.out --error=logs/concat_sort.${pop}.${chr}.%A.err \
    10b_master_gwas_concat_sort.sh ${pop} ${chr}
done
done

for chr in 1 2 3 5 6 8 9 10 11 12 17 19 22 X ; do
for pop in 'metal' ; do
sbatch -p mrivas --qos=high_p --nodes=1 --mem=256000 --cores=6 --time=1-0:00:00 \
    --job-name=concat_sort.${pop}.${chr} --output=logs/concat_sort.${pop}.${chr}.%A.out --error=logs/concat_sort.${pop}.${chr}.%A.err \
    10b_master_gwas_concat_sort.sh ${pop} ${chr}
done
done

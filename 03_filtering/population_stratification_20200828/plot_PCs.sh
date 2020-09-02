#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)

sqc_d="/oak/stanford/groups/mrivas/ukbb24983/sqc"
sqc_v="20200828"
plot_src="$(dirname ${SRCDIR})/pca/plot_PCs.R"


find ${sqc_d}/population_stratification_w24983_${sqc_v} -mindepth 2 -maxdepth 2 -name "*.eigenvec" | sort \
| while read eigenvec_f ; do
for x in $(seq 39) ; do

# test case
# eigenvec_f='/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200828/pca2_covars/ukb24983_s_asian.eigenvec'
# x=1

y=$(expr $x + 1)
out_img=${eigenvec_f}.PC${x}.PC${y}.png

if [ ! -s ${out_img} ] ; then
    echo ${out_img}
    Rscript ${plot_src} ${eigenvec_f} ${out_img} PC${x} PC${y}
fi
done
done

# tar archive

find ${sqc_d}/population_stratification_w24983_${sqc_v} -mindepth 2 -maxdepth 2 -type f -name "*.png" | grep -v __delete  | grep -v .ipynb_checkpoints | sort -V | tar -czvf ${sqc_d}/population_stratification_w24983_${sqc_v}/PCA.plots.tar.gz -T /dev/stdin

readlink -f ${sqc_d}/population_stratification_w24983_${sqc_v}/PCA.plots.tar.gz

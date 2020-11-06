#!/bin/bash
set -beEuo pipefail

BIN='3a_check_results.BINs.20201102-171020.tsv'
BIN_WB='3a_check_results.BINs.WB.20201102-171101.tsv'

cat ${BIN} | egrep -v '#' | awk '($3<100){print $1, $2}' \
| egrep -v 'HC1563|HC1582|HC1584' \
| while read pop GBE_ID ; do 

out_d=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/${pop}-batch
    batch_idxs=$(cat ${BIN_WB} | awk -v pop=${pop} -v GBE_ID=${GBE_ID} '($1==pop && $2 == GBE_ID){print $3}' \
    | sort | comm -23 <(seq 100 | sort) /dev/stdin | tr '\n' ',' | rev | cut -c2- | rev)

    echo ${GBE_ID} ${pop} ${batch_idxs}

    sbatch -p mrivas,normal,owners --nodes=1 \
    --mem=8000 --cores=2 --time=6:00:00 \
    --job-name=gwas.${GBE_ID} \
    --array=${batch_idxs} \
    --output=logs_scratch/gwas.${GBE_ID}.%A.out \
    --error=logs_scratch/gwas.${GBE_ID}.%A.err \
    1_plink.gwas.wrapper.array.sh ${pop} ${GBE_ID}

done
exit 0
#######################################




| egrep -v 'BIN_FC170020544|BIN_FC40001428|HC1018|HC1040|HC1045|HC1061|HC1062|HC1063|HC1064|HC1066|HC1068|HC1069|HC1088|HC1089|HC1090|HC1091|HC1092|HC1093|HC1094|HC1099|HC1101|HC1105|HC1106|HC1112|HC1128|HC1129|HC1130|HC1132|HC1178|HC1179|HC1181|HC1182|HC1183|HC1184|HC1188|HC1192|HC1193|HC1197|HC1198|HC1199|HC1201|HC1202|HC1205|HC1215|HC1216|HC1217|HC1218|HC1224|HC1229|HC1230|HC1234|HC1235|HC1256|HC1257|HC1259|HC1260|HC1265|HC1268|HC1269|HC1271|HC1272|HC1276|HC1277|HC1278|HC1279|HC1283|HC249|HC299|HC305|HC341|HC353|HC369|HC370|HC402|HC740|HC805|HC932|HC938|HC939' \

#!.bin/bash
set -beEu

test_files=(
/oak/stanford/groups/mrivas/dev-ukbb-tools/gwas/ukb_20190327/ukb24983_v2.INI30640.genotyped.glm.linear
/oak/stanford/groups/mrivas/dev-ukbb-tools/gwas/ukb_20190327/ukb24983_v2.BIN23061.genotyped.glm.logistic.hybrid
/oak/stanford/groups/mrivas/dev-ukbb-tools/gwas/ukb_20181109/ukb24983_v2.INI26413.genotyped.glm.linear
/oak/stanford/groups/mrivas/dev-ukbb-tools/gwas/ukb_20181109/ukb24983_v2.BIN21068.genotyped.glm.logistic.hybrid
)

for f in ${test_files[@]} ; do 
    cat $f | python flipfix_and_liftOver.py | head -n5
    echo ""
    cat $f | python flipfix_and_liftOver.py --to38 | head -n5
    echo ""
done


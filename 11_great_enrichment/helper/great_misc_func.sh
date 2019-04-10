get_onto_list () {
    assembly=$1
    onto_file="$(dirname $(dirname $(readlink -f $0)))/misc/great.ontology.${assembly}.20171029.txt"

    cat ${onto_file} | egrep -v '^#' | cut -f1 | tr "\n" "," | rev | cut -c2- | rev
}

great_cols () {
    echo "ID Desc DAGLvl BRank BPval BBonf BFDR BFold BExp Bn Bk BProb RegionSetCov HRank HPval HBonf HFDR HFold HExp HN HK Hn Hk GeneSetCov TermGeneCov Genes Regions" | tr " " "\t"
}

cat_or_zcat () {
    file=$1
    if [ ${file%.gz}.gz == ${file} ] ; then zcat ${file} ; else cat ${file} ; fi
}

get_col_by_key () {
    file=$1
    key=$2
    cat_or_zcat $file | awk 'NR==1' | tr "\t" "\n" | awk -v key=$key '($0==key){print NR}'
}

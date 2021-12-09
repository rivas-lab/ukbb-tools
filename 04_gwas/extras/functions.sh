cat_or_zcat () {
    local file=$1
    if [ "${file%.gz}.gz" == "${file}" ] || [ "${file%.bgz}.bgz" == "${file}" ] ; then
        zcat "${file}"
    elif [ "${file%.zst}.zst" == "${file}" ] ; then
        zstdcat "${file}"
    else
        cat "${file}"
    fi
}


find_compressed_file_if_exist () {
    local file=$1
    if [ -s "${file}.zst" ] ; then
        echo "${file}.zst"
    elif [ -s "${file}.gz" ] ; then
        echo "${file}.gz"
    else
        echo "${file}"
    fi
}


cat_or_zcat_prefix () {
    local file=$1
    if  [ -s "${file}.zst" ] ||
        [ -s "${file}.gz" ] ||
        [ -s "${file}" ] ; then
        cat_or_zcat $(find_compressed_file_if_exist "${file}")
    fi
}


load_plink2 () {
    # We installed PLINK2 software as a software module in our HPC system.
    # We call `ml load plink2` for the specified version. By doing this,
    # it updates the PATHs so that we can execute plink2 software.

    plink2_version=$1

    check_avx2_flag=$(cat /proc/cpuinfo | grep flags | grep -i avx2 | uniq | cat /dev/stdin <(echo avx2) | grep -i avx2 | wc -l)

    if [ "${check_avx2_flag}" -eq 2 ] ; then
        ml load "plink2/${plink2_version}"
    else
        ml load "plink2/${plink2_version}-non-AVX2"
    fi
}


get_phe_file () {
    local GBE_ID=$1
    cat /oak/stanford/groups/mrivas/users/${USER}/repos/rivas-lab/ukbb-tools/05_gbe/array-combined/phenotype_info.tsv \
    | awk -v GBE_ID="${GBE_ID}" '($1 == GBE_ID){print $NF}'
}


get_plink_suffix () {
    local GBE_ID=$1

    GBE_CAT=$(echo $GBE_ID | sed -e "s/[0-9]//g")

    if [ "${GBE_CAT}" == "QT_FC" ] || [ "${GBE_CAT}" == "INI" ] ; then
        echo "glm.linear"
    else
        echo "glm.logistic.hybrid"
    fi
}


get_covar_PCs () {
    local pop=$1

    if [ "${pop}" == "others" ] || [ "${pop}" == "related" ] ; then
        echo "Global_PC1-Global_PC18"
    else
        echo "PC1-PC10"
    fi
}


show_scracth_if_exists () {
    local f=$1
    scratch_f=$(echo $f | sed -e "s%/oak/stanford/%/scratch/%g")
    if [ -f "${scratch_f}" ] ; then echo "${scratch_f}" ; else echo "${f}" ; fi
}


############################################################
# functions to combine files
############################################################


combine_get_file_name_with_batch_idx () {
    template_f=$1
    batch_idx=$2
    echo ${template_f} | sed -e "s/__BATCH__/${batch_idx}/g"
}


combine_check_files () {
    template_f=$1
    n_batch=$2

    seq "${n_batch}" | tr ' ' '\n' | while read -r batch_idx ; do
        f=$(combine_get_file_name_with_batch_idx ${template_f} ${batch_idx})
        if [ ! -s "${f}" ] && [ ! -s "${f}.zst" ] && [ ! -s "${f}.gz" ] ; then
            echo "[missing] ${batch_idx} ${f}"
        fi
    done
}


combine_get_and_cat_file_with_batch_idx () {
    template_f=$1
    batch_idx=$2

    f=$(combine_get_file_name_with_batch_idx ${template_f} ${batch_idx})

    if  [ -s "${f}.zst" ] ||
        [ -s "${f}.gz" ] ||
        [ -s "${f}" ] ; then
        cat_or_zcat $(find_compressed_file_if_exist "${f}")
    fi
}


combine_log_files () {
    # combine the log files
    log_template_f=$1
    n_batch=$2
    log_combined_f=$3
    cores=1
    if [ $# -gt 3 ] ; then cores=$4 ; fi

    if [ ! -d $(dirname ${log_combined_f}) ] ; then mkdir -p $(dirname ${log_combined_f}) ; fi

    seq ${n_batch} | tr ' ' '\n' | while read -r batch_idx ; do
        if [ "${batch_idx}" -ne 1 ] ; then echo "" ; fi
        f=$(combine_get_file_name_with_batch_idx ${log_template_f} ${batch_idx})
        echo "## ${f}"
        cat_or_zcat ${f}
    done | bgzip -l9 -@${cores} > ${log_combined_f%.gz}.gz
}


combine_plink_files () {
    # combine the log files
    plink_template_f=$1
    n_batch=$2
    plink_combined_f=$3
    cores=1
    if [ $# -gt 3 ] ; then cores=$4 ; fi

    if [ ! -d $(dirname ${plink_combined_f}) ] ; then mkdir -p $(dirname ${plink_combined_f}) ; fi

    {
        combine_get_and_cat_file_with_batch_idx ${plink_template_f} 1 | awk 'NR==1'

        seq ${n_batch} | tr ' ' '\n' | while read -r batch_idx ; do
            combine_get_and_cat_file_with_batch_idx ${plink_template_f} ${batch_idx} | awk 'NR>1'
#         done | sort -k1,1V -k2,2n -k3,3 # the file should be sorted..
        done
    } | bgzip -l9 -@${cores} > ${plink_combined_f%.gz}.gz
}

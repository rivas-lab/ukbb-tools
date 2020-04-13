#!/bin/bash
set -beEuo pipefail

usage () {
    echo "$0: Updates PLINK with the version specified as the first argument."
    echo "usage: bash $0 MMMMYYDD"
}

if [ $# -lt 1 ] ; then usage >&2 ; exit 1 ; fi

plink_version=$1

get_plink_url () {
    local os=$1
    local version=$2

    local PLINK_URL="http://s3.amazonaws.com/plink2-assets/plink2_{{ OS }}_{{ VERSION }}.zip"

    echo ${PLINK_URL} | sed -e "s/{{ OS }}/${os}/g" | sed -e "s/{{ VERSION }}/${version}/g"
}

install_plink () {
    local os=$1
    local version=$2
    if [ $# -gt 2 ] ; then set_default=$3 ; else set_default=1 ; fi

    local install_dir_root="/oak/stanford/groups/mrivas/software/plink2"
    local module_dir_root="/home/groups/mrivas/.modules/plink2"

    if [ ${os} == "linux_x86_64" ] ; then local suffix="-non-AVX2" ; else local suffix="" ; fi
    local install_d=${install_dir_root}/${plink_version}${suffix}
    local plink_url=$(get_plink_url ${os} $plink_version)
    local module_f=${module_dir_root}/${plink_version}${suffix}.lua
    local module_template=${module_dir_root}/template-${os}.txt

    if [ ! -f ${module_template} ] ; then echo "${os} is not supported in this auto install script" >&2 ; exit 1; fi

    mkdir -p ${install_d} && cd ${install_d}
    wget ${plink_url}
    unzip $(basename ${plink_url})
    chmod 770 plink2
    cd -
    cat ${module_template} | sed -e "s/{{ VERSION }}/${version}/g" > ${module_f}
    if [ ${set_default} -eq 0 ] ; then
        cd ${module_dir_root}
        ln -sf ${module_f} default
        cd -
    fi
}

install_plink linux_x86_64 $plink_version 1
install_plink linux_avx2   $plink_version 0


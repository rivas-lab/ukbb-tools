#!/bin/bash
set -beEuo pipefail

version=$1

repo_url="git@github.com:rivas-lab/ukbb-tools.git"

software_dir="/oak/stanford/groups/mrivas/software/ukbb-tools"

module_dir="/home/groups/mrivas/.modules/ukbb-tools"


cd ${software_dir} 

git clone ${repo_url}

cd ukbb-tools
git tag -a ${version} -m "freeze ${version}"
git push origin master --tags
cd ..

mv ukbb-tools ${version}

cat ${module_dir}/module.template.txt | sed -e "s/__VERSION__/${version}/g" > ${module_dir}/${version}.lua

cd ${module_dir}
ln -sf ${version}.lua default


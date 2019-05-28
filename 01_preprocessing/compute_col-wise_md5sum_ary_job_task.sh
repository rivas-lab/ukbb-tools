#!/bin/bash
set -beEuo pipefail

script="$(dirname $(readlink -f $0))/compute_col-wise_md5sum.sh"

in_file=$(readlink -f $1)
task_id=$2
batch_size=$3
shift 3

(( end_idx = task_id * batch_size ))
(( start_idx = end_idx - batch_size + 1 ))

bash $script -s ${start_idx} -e ${end_idx} ${in_file} $@


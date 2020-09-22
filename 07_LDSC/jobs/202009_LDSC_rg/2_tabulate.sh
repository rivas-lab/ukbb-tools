#!/bin/bash
set -beEuo pipefail

find /oak/stanford/groups/mrivas/ukbb24983/array-combined/ldsc/rg/ -type f -name "*.log" \
    | bash /oak/stanford/groups/mrivas/users/${USER}/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ldsc_rg_view.sh -l /dev/stdin


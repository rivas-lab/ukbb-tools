#!/bin/bash
set -beEuo pipefail

{
    echo "#wc_l file"
    find data -name "*glm*" | sort | while read f ; do echo $(cat $f | wc -l) $f ; done
} | tr ' ' '\t' > 4_wc_l.tsv

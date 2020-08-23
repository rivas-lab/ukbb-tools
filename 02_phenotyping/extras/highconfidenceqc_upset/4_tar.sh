#!/bin/bash
set -beEuo pipefail

data_d='/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/highconfidenceqc_upset/upset_plot'

find ${data_d} -name "*.png" | tar --transform 's/.*\///g' --group root --owner root -czvf ${data_d}/upset_plot.tar.gz -T /dev/stdin

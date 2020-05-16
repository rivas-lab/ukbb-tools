#!/bin/bash

cat original_hc_gbe_map.tsv <(awk -F'\t' '{if (NR > 1) {print $9"\t"$3}}' TTE_HC.tsv) <(awk -F'\t' '{if (NR > 1) {print $9"\t"$3}}' AD_HC.tsv) > highconfidenceqc_gbe_map.tsv

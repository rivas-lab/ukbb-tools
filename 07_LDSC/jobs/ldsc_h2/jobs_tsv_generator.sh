#!/bin/bash
set -beEuo pipefail

sumstats_root=$OAK/dev-ukbb-tools/gwas

find $sumstats_root -mindepth 2 -maxdepth 2 -type f -name "ukb24983_v2.*.glm.*" \
	| egrep -v "id$" \
	| sort > jobs.tsv


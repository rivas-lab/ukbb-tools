#!/bin/bash
set -beEuo pipefail

cat /oak/stanford/groups/mrivas/users/${USER}/repos/rivas-lab/ukbb-tools/05_gbe/phenotype_info.tsv \
| cut -f1 | sort -V > GBE_ID.lst


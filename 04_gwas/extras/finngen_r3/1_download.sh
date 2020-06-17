#!/bin/bash
set -beEuo pipefail

cd /scratch/groups/mrivas/public_data/summary_stats/finngen_r3

gsutil -m cp -r gs://finngen-public-data-r3/summary_stats/ .
gsutil -m cp -r gs://finngen-public-data-ld/imputation_panel_v1/ .
gsutil -m cp -r gs://finngen-public-data-r3/finemapping/ .

#!/bin/bash
set -beEuo pipefail

sbatch -p mrivas --qos=high_p --nodes=1 --mem=8000 --cores=1 --time=7-0:00:00 \
    --job-name=dl --output=logs/dl.%A_%a.out --error=logs/dl.%A_%a.err \
    --array=2,6,15,20,21 download.sh


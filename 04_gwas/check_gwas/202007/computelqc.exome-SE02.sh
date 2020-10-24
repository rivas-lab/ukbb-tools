#!/bin/bash
set -beEuo pipefail

ml load R/3.6 gcc

Rscript computelqc.exome-SE02.R $@

#!/bin/bash
set -beEuo pipefail

ml load plink2

d="/oak/stanford/groups/mrivas/ukbb24983"

plink2 --memory 60000 --threads 6 \
--pfile ${d}/array_combined/pgen/ukb24983_cal_hla_cnv vzs \
--extract /oak/stanford/groups/mrivas/ukbb24983/cal/coding-only/coding.vars.lst \
--make-pgen vzs \
--out ${d}/cal/coding-only/ukb24983_cal_coding

mv ${d}/cal/coding-only/ukb24983_cal_coding.log ${d}/cal/coding-only/ukb24983_cal_coding.pgen.log

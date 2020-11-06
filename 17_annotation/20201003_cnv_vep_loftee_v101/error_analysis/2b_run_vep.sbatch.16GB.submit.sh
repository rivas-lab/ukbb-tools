#!/bin/bash
set -beEuo pipefail

log_d=logs

if [ ! -d ${log_d} ] ; then mkdir -p ${log_d} ; fi

ml load resbatch

# sbatch -p mrivas,normal,owners \
sbatch -p mrivas,normal \
--time=2-0:0:00 --mem=16000 --nodes=1 --cores=1 --job-name=vep --output=${log_d}/vep.%A_%a.out --error=${log_d}/vep.%A_%a.err \
--array=2,3,6,7,8,10,21,29,37,38,39,42,44,51,56,60,67,70,72,76,89,112,114,115,126,128,130,134,138,160,165,177,181,185,187,188,190,193,199,225,230,231,234,237,246,247,256,258,264,279,280,281,282,284,287,288,296,301,304,307,313,321,344,346,351,352,359,365,368,382,385,390,407,408,409,414,418,420,421,425,429,432,433,437,442,444,454,455,457,473,476,479,484,486,489,491,500,504,508,509,512,513,514,521,523,556,557,558,572,581,590,591,602,609,612,616,620,621,631,644,648,651,654,667,668,672,676,684,695,698,701,704,717,724,728,729,734,743,757,761,762,767,777,782,783,788,789,798,818,819,823,824,832,833,844,846,848,851,852,854,855,857,861,866,867,869,870,872,876,877,906,924,955,956,964,966,971,977,978,980,982,983,985,988,993 \
${parallel_sbatch_sh} 2_run_vep.sh ${parallel_idx} 1

exit 0
/oak/stanford/groups/mrivas/ukbb24983/cnv/pgen/cnv.pvar
# 275181 lines, 275180 variants

998 batches
276 variants each
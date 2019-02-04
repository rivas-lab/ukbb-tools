#!/bin/bash
set -beEu -o pipefail

if [ $# -lt 1 ] ; then echo "usage: $0 <err_dir> [job_num=4]" >&2 ; exit 1; fi

err_dir=$1

if [ $# -lt 2 ] ; then job_num=4 ; else job_num=$2 ; fi

parallel --jobs=$job_num "cat {} | grep 'array-start' | rev | awk '{print \$1}' | rev" ::: $( find $err_dir -name "*.err" ) \
	| sort -g


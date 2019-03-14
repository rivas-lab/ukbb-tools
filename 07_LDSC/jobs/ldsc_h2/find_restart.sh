#!/bin/bash
set -beEu -o pipefail

if [ $# -lt 2 ] ; then echo "usage: $0 <n_tasks> <err_dir> [job_num=4]" >&2 ; exit 1; fi

n_tasks=$1
err_dir=$2

if [ $# -lt 3 ] ; then job_num=4 ; else job_num=$3 ; fi

comm -13 <( bash find_started.sh $err_dir $job_num | sort ) <( seq 1 $n_tasks | sort ) \
	| sort -g | tr "\n" "," | rev | cut -c2- | rev


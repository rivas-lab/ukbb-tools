#!/bin/bash
#set -beEuo pipefail

field_id=23155

echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname=$(hostname); SLURM_JOBID=${SLURM_JOBID:=0}; SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:=24}" >&2

c=${SLURM_ARRAY_TASK_ID:=22}
#./gfetch ${field_id} -c$c -m
./gfetch ${field_id} -c$c

echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname=$(hostname); SLURM_JOBID=${SLURM_JOBID}; SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2


#!/bin/bash
set -beEuo pipefail

in_var_list=$1

one_array_list="/oak/stanford/groups/mrivas/ukbb24983/sqc/one_array_variants.txt"

cat ${in_var_list} | sort | comm -12 /dev/stdin <(cat $one_array_list | sort) > ${in_var_list%.lst}.one_array_variants.lst

cat ${in_var_list} | sort | comm -23 /dev/stdin <(cat $one_array_list | sort) > ${in_var_list%.lst}.both_arrays_variants.lst

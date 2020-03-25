#!/bin/bash
set -beEuo pipefail

seq 1 22 | parallel -j4 --eta -k 'ukbgene hap -c{} -a.ukbkey'
seq 1 22 | parallel -j4 --eta -k 'ukbgene hap -c{} -a.ukbkey -m'


#!/bin/bash
set -beEuo pipefail

cd /oak/stanford/groups/mrivas/users/ytanigaw/sandbox/20200612_metal_power_analysis

for npop in 4 5 6 7 ; do
    cat INI50.lst | awk -v npop=${npop} 'NR <= npop' | 
    bash ~/repos/rivas-lab/ukbb-tools/18_metal/run_metal.sh -o INI50.${npop}pops -f /dev/stdin
done

for npop in 5 6 ; do
    cat INI50.lst | grep -v e_asian | awk -v npop=${npop} 'NR <= npop' | 
    bash ~/repos/rivas-lab/ukbb-tools/18_metal/run_metal.sh -o INI50.${npop}pops-noEA -f /dev/stdin
done

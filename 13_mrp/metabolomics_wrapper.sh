#!/bin/bash

find /oak/stanford/groups/mrivas/users/ytanigaw/data/.20200605_NGH_UKB_phase1_metabolomics/gwas_exome/20201104 -name "*gz" | cut -d'.' -f3 | sort -u  | grep -v NMR_metabolomics > metabolomics_phenos.tsv

for phe in $(cat metabolomics_phenos.tsv); do
    echo $phe
    #sbatch -p owners,mrivas --mem=128000 --nodes=1 --cores=8 -t 2-00:00:00 mrp_rv_metabolomics.sh $phe;
    sbatch -p owners,mrivas --mem=128000 --nodes=1 --cores=8 -t 2-00:00:00 mrp_rv_metabolomics_mpc_pli.sh $phe;
done

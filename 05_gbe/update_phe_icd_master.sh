#/bin/bash
echo "Updating phenotype_info.tsv..."
bash update_phe_info.sh
echo "Updated phenotype_info.tsv."
echo "Updating icdinfo.txt..."
bash update_icdinfo.sh
echo "Updated icdinfo.txt."
echo "Submitting job to generate exome_phenotype_info.tsv..."
sbatch -p normal,mrivas,owners --mem=64000 -t 1-00:00:00 -J exm_phen_info --wrap="python make_exome_phenotype_info.py"
echo "Submitting job to generate master.phe..."
sbatch -p normal,mrivas,owners --mem=64000 -t 1-00:00:00 -J masterPhe --wrap="python combine_phe.py"
echo "master.phe is brewing! Please check back in 1 business day."

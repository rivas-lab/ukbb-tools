#/bin/bash

bash update_phe_info.sh
bash update_gbe_icdinfo.sh
sbatch -p normal,mrivas,owners --mem=64000 -t 1-00:00:00 -J masterPhe --wrap="python combine_phe.py"

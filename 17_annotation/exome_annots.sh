#!/bin/bash

sbatch --mem=64000 -t 3-00:00:00 -p mrivas -J sas --nodes=1 --cores=8 annotate_bims.sh ukb_exm_spb-s_asian.bim GRCh38
sbatch --mem=64000 -t 3-00:00:00 -p mrivas -J eas --nodes=1 --cores=8 annotate_bims.sh ukb_exm_spb-e_asian.bim GRCh38
sbatch --mem=64000 -t 3-00:00:00 -p mrivas -J wb --nodes=1 --cores=8 annotate_bims.sh ukb_exm_spb-white_british.bim GRCh38
sbatch --mem=64000 -t 3-00:00:00 -p mrivas -J nbw --nodes=1 --cores=8 annotate_bims.sh ukb_exm_spb-non_british_white.bim GRCh38
sbatch --mem=64000 -t 3-00:00:00 -p mrivas -J african --nodes=1 --cores=8 annotate_bims.sh ukb_exm_spb-african.bim GRCh38

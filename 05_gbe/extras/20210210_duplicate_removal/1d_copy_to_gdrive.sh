#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})

source 0_parameters.sh

ml load rclone

rclone copy ${potential_dups_GBE_w_md5_f} gdrive://${gdrive_path}
exit 0



rclone copy ${phenotyping_annotation_f} gdrive://${gdrive_path}

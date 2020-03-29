#!/bin/bash
set -beEuo pipefail

foo () {

#     local plink_common_opts=" --threads ${cpu} --memory ${mem} "
    local bgen_dir="/scratch/groups/mrivas/ukbb24983/hap/download"
#     local sample="${bgen_dir}/ukb24983_hap_chr${c}_v2_s${s_cnt}.sample"    
#     local scratch_dir="/scratch/groups/mrivas/ukbb24983/hap/pgen"
#     local oak_dir="/oak/stanford/groups/mrivas/ukbb/24983/hap/pgen"
#     local base="ukb_hap_chr${c}_v2"
#     local oak_f="${oak_dir}/${base}"
#     local scr_f="${scratch_dir}/${base}"
#     local bgen=${bgen_dir}/${base}.bgen
#     local ref_fa="/scratch/groups/mrivas/public_data/genomes/hg19/hg19.fa"

    cat-bgen -clobber -g $(for c in $(seq 1 21) ; do echo "${bgen_dir}/ukb_hap_chr${c}_v2.bgen" ; done) -og ${bgen_dir}/ukb_hap_chrAUTO_v2.bgen
}

foo

# [ytanigaw@sh02-09n53 ~/repos/rivas-lab/ukbb-tools/03_filtering/hap]$ bash 3_concat.sh

# Welcome to cat-bgen
# (version: 44fcabbc5c3892e2241fca593db789e0d97e92ed)

# (C) 2009-2017 University of Oxford

# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr1_v2.bgen" (1 of 21, 53260 variants)...
# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr2_v2.bgen" (2 of 21, 52631 variants)...
# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr3_v2.bgen" (3 of 21, 44263 variants)...
# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr4_v2.bgen" (4 of 21, 41213 variants)...
# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr5_v2.bgen" (5 of 21, 38734 variants)...
# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr6_v2.bgen" (6 of 21, 45886 variants)...
# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr7_v2.bgen" (7 of 21, 36061 variants)...
# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr8_v2.bgen" (8 of 21, 33741 variants)...
# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr9_v2.bgen" (9 of 21, 29133 variants)...
# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr10_v2.bgen" (10 of 21, 32626 variants)...
# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr11_v2.bgen" (11 of 21, 33134 variants)...
# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr12_v2.bgen" (12 of 21, 31429 variants)...
# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr13_v2.bgen" (13 of 21, 22134 variants)...
# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr14_v2.bgen" (14 of 21, 21337 variants)...
# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr15_v2.bgen" (15 of 21, 20859 variants)...
# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr16_v2.bgen" (16 of 21, 23774 variants)...
# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr17_v2.bgen" (17 of 21, 22215 variants)...
# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr18_v2.bgen" (18 of 21, 19250 variants)...
# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr19_v2.bgen" (19 of 21, 19139 variants)...
# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr20_v2.bgen" (20 of 21, 17197 variants)...
# Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr21_v2.bgen" (21 of 21, 9793 variants)...
# Finished writing "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chrAUTO_v2.bgen" (487409 samples, 647809 variants).

# Thank you for using cat-bgen.

[ytanigaw@sh02-09n53 ~/repos/rivas-lab/ukbb-tools/03_filtering/hap]$ bash 3_concat.sh

Welcome to cat-bgen
(version: 44fcabbc5c3892e2241fca593db789e0d97e92ed)

(C) 2009-2017 University of Oxford

Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr1_v2.bgen" (1 of 22, 53260 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr2_v2.bgen" (2 of 22, 52631 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr3_v2.bgen" (3 of 22, 44263 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr4_v2.bgen" (4 of 22, 41213 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr5_v2.bgen" (5 of 22, 38734 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr6_v2.bgen" (6 of 22, 45886 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr7_v2.bgen" (7 of 22, 36061 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr8_v2.bgen" (8 of 22, 33741 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr9_v2.bgen" (9 of 22, 29133 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr10_v2.bgen" (10 of 22, 32626 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr11_v2.bgen" (11 of 22, 33134 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr12_v2.bgen" (12 of 22, 31429 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr13_v2.bgen" (13 of 22, 22134 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr14_v2.bgen" (14 of 22, 21337 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr15_v2.bgen" (15 of 22, 20859 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr16_v2.bgen" (16 of 22, 23774 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr17_v2.bgen" (17 of 22, 22215 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr18_v2.bgen" (18 of 22, 19250 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr19_v2.bgen" (19 of 22, 19139 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr20_v2.bgen" (20 of 22, 17197 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr21_v2.bgen" (21 of 22, 9793 variants)...
Adding file "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr22_v2.bgen" (22 of 22, 10911 variants)...
Error: input file #22 ( "/scratch/groups/mrivas/ukbb24983/hap/download/ukb_hap_chr22_v2.bgen" ) has the wrong flags (80000009, expected 9).  Quitting.

Thank you for using cat-bgen.
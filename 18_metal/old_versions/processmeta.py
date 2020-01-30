from __future__ import print_function 

import os
import sys


trait = sys.argv[1]

fout = open(trait + '.masterfile','w')

print("# === DESCRIBE AND PROCESS THE INPUT INPUT FILES ===", file = fout)
print("SCHEME STDERR", file = fout)
print("MARKERLABEL ID", file = fout)
print("ALLELE A1 REF",file = fout)
print("EFFECT BETA", file = fout)
print("PVALUE P", file = fout)
print("STDERR SE", file = fout)
print("WEIGHT OBS_CT", file = fout)
print("\n", file = fout)
print("\n", file = fout)
#e_asian
#non_british_white
#s_asian
#african
#w_british
#ukb24983_v2_hg19.Apolipoprotein_B_adjust_statins.genotyped.glm.linear.gz
print("PROCESS /oak/stanford/groups/mrivas/projects/biomarkers_rivas/main/w_british/ukb24983_v2_hg19." + trait + ".genotyped.glm.linear.gz", file = fout)
print("PROCESS /oak/stanford/groups/mrivas/projects/biomarkers_rivas/main/non_british_white/ukb24983_v2_hg19." + trait + ".genotyped.glm.linear.gz", file = fout)
print("PROCESS /oak/stanford/groups/mrivas/projects/biomarkers_rivas/main/s_asian/ukb24983_v2_hg19." + trait + ".genotyped.glm.linear.gz", file = fout)
print("PROCESS /oak/stanford/groups/mrivas/projects/biomarkers_rivas/main/african/ukb24983_v2_hg19." + trait + ".genotyped.glm.linear.gz", file = fout)
# Excluded East Asian
print("\n", file = fout)
print("OUTFILE METAANALYSIS_" + trait + "_ .tbl", file = fout)
print("MINWEIGHT 10000", file = fout)
print("\n", file = fout)
#print("ANALYZE", file = fout)
print("ANALYZE HETEROGENEITY", file = fout)
#print("SCHEME   STDERR", file = fout)
print("\n", file = fout)
print("QUIT", file = fout)


##CHROM  POS     ID      REF     ALT     A1      TEST    OBS_CT  BETA    SE      T_STAT  P

#names ID      SNP     TOTAL_MAF       CHROM   POS     REF     ALT1    ALT     A1      AX      A1_CT   ALLELE_CT       A1_FREQ MACH_R2 TEST    OBS_CT  BETA    SE
#      T_STAT  P
#


fout.close()




def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)
    

job_directory = "%s/.job" %os.getcwd()

# Make top level directories
mkdir_p(job_directory)

job_file = os.path.join(job_directory,"%s.job" %trait)


with open(job_file,'w') as fh:
    fh.writelines("#!/bin/bash\n")
    fh.writelines("#SBATCH --job-name=%s.job\n" % trait)
    fh.writelines("#SBATCH --output=.job/%s.out\n" % trait)
    fh.writelines("#SBATCH --error=.job/%s.err\n" % trait)
    fh.writelines("#SBATCH --time=1-00:00\n")
    fh.writelines("#SBATCH --mem=48000\n")
    fh.writelines("#SBATCH -p owners,normal,mrivas\n")
    fh.writelines("/scratch/PI/mrivas/users/hanna/bin/generic-metal/metal %s.masterfile \n" %trait)
    fh.close()

os.system("sbatch %s" %job_file)

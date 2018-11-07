#!/bin/python

_README_='''
A script for running GWAS with UK Biobank data (array and/or imputed genotypes) using PLINK.

Author: Matthew Aguirre (SUNET: magu)
'''

def make_plink_command(bpFile, pheFile, outFile, pop, related=False, plink1=False, arrayCovar=False):
    # paths to plink genotypes, input phenotypes, output directory are passed
    # paths to qc files (ethnic groups, related individual list) are coded here
    covarFile='/oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar.phe'
    return " ".join(["plink" if plink1 else "plink2", 
                     "--bfile" if plink1 else "--pfile", bpFile,
                     "--pheno", pheFile, "--pheno-quantile-normalize",
                     "--glm firth-fallback hide-covar",
                     "--keep /path/to/ethnic/group",
                     "--exclude /latest/redacted.list /path/to/related/inds.list",
                     "--extract {0}".format("on_both.arrays" if not arrayCovar else "on_one.array"),
                     "--covar", covarFile, 
                     "--covar-name age sex {0} PC1-PC4".format("Array" if arrayCovar else ""), 
                     "--out {0}".format(outFile)]) 

def make_batch_file(batchFile, plinkCmd, memory, time, partitions):
    with open(batchFile, 'w') as f:
       f.write("\n".join(["#!/bin/bash","",
                          "#SBATCH --job-name=RL_GWAS",
                          "#SBATCH --output={}".format(os.path.join(os.path.dirname(batchFile), "rl-gwas_%A-%a.out")),
                          "#SBATCH --memory={}".format(memory),
                          "#SBATCH --time={}".format(time),
                          "#SBATCH -p {}".format(','.join(partitions)),
                          "", plinkCmd]))
    return batchFile


def run_gwas(kind, pheFile, outDir='', pop='white_british', related=False, plink1=False, 
             logDir='', memory="24000", time="1-00:00:00", partition=["normal","owners"]):
    import os
    # ensure usage
    if not os.path.isfile(pheFile):
        raise ValueError("Error: phenotype file {0} does not exist!".format(pheFile))
    if not os.path.isdir(outDir):
        raise ValueError("output directory {0} does not exist!".format(outDir))
    if not os.path.isdir(logDir):
        print("Warning: Logging directory either does not exist or was unspecified. Defaulting to {0}...".format(os.getcwd()))
        logDir=os.getcwd()
    if pop not in ['all', 'white_british', 'african', 's_asian', 'e_asian']:
        raise ValueError("population must be one of (all, white_british, african, s_asian, e_asian)")
    # derive for the below
    imp_bfile_path='imp.pgen'
    cal_bfile_path='cal.pgen'
    outFile='out'
    if kind == 'imputed':
        # make this an array job with nJobs=22
        cmd = make_plink_command(bpFile  = imp_bfile_path,
                                 pheFile = pheFile,
                                 outFile = outFile,
                                 pop     = pop,
                                 related = related,
                                 plink1  = plink1,
                                 arrayCovar = True)
    elif kind == 'genotyped': 
        # need a run with array a covariate, and one without
        cmd1 = make_plink_command(bpFile  = cal_bfile_path,
                                  pheFile = pheFile,
                                  outFile = outFile,
                                  pop     = pop,
                                  related = related,
                                  plink1  = plink1,
                                  arrayCovar = True)
        cmd2 = make_plink_command(bpFile  = cal_bfile_path,
                                  pheFile = pheFile,
                                  outFile = outFile,
                                  pop     = pop,
                                  related = related,
                                  plink1  = plink1,
                                  arrayCovar = False)
        cmd = "\n".join([cmd1, cmd2, # this is the plink part, below joins the two files
                         "cat out1 out2 | sort -k1,1n -kN,lookItUp > finalout",
                         "cat out1.log out2.log > finalout.log",
                         "rm out1.* out2.*" # this is problematic if one is a substring, so program 1/2 yourself to avoid
                        ])
    else:
        raise ValueError("argument kind must be one of (imputed, genotyped): {0} was provided".format(kind))
    sbatch = make_batch_file(batchFile = os.path.join(logDir, "some_name_{}.sh".format(kind)),
                             plinkCmd  = cmd,
                             memory    = memory,
                             time      = time,
                             partitions = partition)
    return
    os.system(" ".join("sbatch", 
                       "{0} -o {1}".format("--array=1-22" if imputed else "", outDir),
                       sbatch))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=_README_
    )
    parser.add_argument('--plink1', dest='plink1', action='store_true',
                            help='Flag to run GWAS with PLINK v1.9 instead of v2.0')
    parser.add_argument('--run-array', dest="arr", action='store_true',
                            help='Run GWAS on directly genotyped (array) data')
    parser.add_argument('--run-imputed', dest="imp", action='store_true',
                            help='Run GWAS on imputed data') 
    parser.add_argument('--pheno', dest="pheno", required=True, nargs='*',
                            help='Path to phenotype file(s)')
    parser.add_argument('--out', dest="outDir", required=True, nargs=1,
                            help='Path to desired output *directory*. Summary stats will be output according to phenotype name (derived from passed file) and Rivas Lab specification for GBE. Defaults to current working directory.')
    parser.add_argument('--population', dest="pop", required=False, default="white_british", nargs=1,
                            help='Flag to indicate which ethnic group to use for GWAS. Must be one of all, white_british, e_asian, s_asian, african')
    parser.add_argument('--keep-related', dest="relatives", action='store_true',
                            help='Flag to keep related individuals in GWAS. Default is to remove them.')
    parser.add_argument('--batch-memory', dest="sb_mem", required=False, default="16000", nargs=1,
                            help='For underlying batch job submission: Amount of memory (in MB) to request. Default is 16000.')
    parser.add_argument('--batch-time', dest="sb_time", required=False, default="24:00:00", nargs=1,
                            help='For underlying batch job submission: Amount of time (DD-HH:MM:SS) to request. Default is 24 hours.')
    parser.add_argument('--batch-partitions', dest="sb_parti", required=False, default=["normal","owners"], nargs='*',
                            help='For underlying batch job submission: Compute partition to submit jobs. Default is normal,owners.')
    parser.add_argument('--log-dir', dest="log", required=False, nargs=1,
                            help="Directory in which to place log files from this script and its related SLURM jobs. Default is current working directory.")
    args = parser.parse_args()
    # feature add: genotype model  
    print(args) 
    import os
    # ensure usage:
    if not args.arr and not args.imp:
        raise ValueError("Error: at least one of --run-array, --run-imputed must be passed")
    # lol i hope this works
    if args.imp:
        run_gwas(kind='imputed', pheFile=args.pheno[0], outDir=args.outDir[0], pop=args.pop[0], related=args.relatives, plink1=args.plink1, 
                 logDir=args.log[0], memory=args.sb_mem[0], time=args.sb_time[0], partition=args.sb_parti)
    if args.arr:
        run_gwas(kind='genotyped', pheFile=args.pheno[0], outDir=args.outDir[0], pop=args.pop[0], related=args.relatives, plink1=args.plink1, 
                 logDir=args.log[0], memory=args.sb_mem[0], time=args.sb_time[0], partition=args.sb_parti)
    
 

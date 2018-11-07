#!/bin/python

_README_='''
A script for running GWAS with UK Biobank data (array and/or imputed genotypes) using PLINK.

Author: Matthew Aguirre (SUNET: magu)
'''

def make_plink_command(bpFile, pheFile, outFile, pop, related=False, plink1=False, arrayCovar=False):
    # paths to plink genotypes, input phenotypes, output directory are passed
    # paths to qc files (covars, ethnic groups, related individuals, ) are coded here
    covarFile='/oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar.phe'
    # make sure this includes most recent update of redacted individuals!
    popFile='/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_{0}.phe'.format(pop) if pop != 'all' else ''
    unrelatedFile='/oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_v2.used_in_pca.phe' if not related else ''
    # infer whether we need to subset variant lists to have array as a covariate 
    # (only true for genotyped runs)
    if arrayCovar and '/cal/' in bpFile: 
        arrayVarFile='/oak/stanford/groups/mrivas/ukbb24983/sqc/both_array_variants.txt'
    if not arrayCovar and '/cal/' in bpFile:
        arrayVarFile='/oak/stanford/groups/mrivas/ukbb24983/sqc/one_array_variants.txt'
    else:
        arrayVarFile=''
    # paste together the command from constituent parts
    return " ".join(["plink" if plink1 else "plink2", 
                     "--bfile", bpFile, "--chr 1-22",
                     "--pheno", pheFile, "--pheno-quantile-normalize",
                     "--glm firth-fallback hide-covar",
                     "--keep {0}".format(popFile) if popfile else "", 
                     "--remove {0}".format(unrelatedFile) if unrelatedFile else "",
                     "--extract {0}".format(arrayVarFile) if arrayVarFile else "",
                     "--covar", covarFile, 
                     "--covar-name age sex", "Array" if arrayCovar else "", "PC1-PC4",
                     "--out", outFile]) 

def make_batch_file(batchFile, plinkCmd, memory, time, partitions):
    with open(batchFile, 'w') as f:
       f.write("\n".join(["#!/bin/bash","",
                          "#SBATCH --job-name=RL_GWAS",
                          "#SBATCH --output={}".format(os.path.join(os.path.dirname(batchFile), "rl-gwas.%A-%a.out")),
                          "#SBATCH --mem={}".format(memory),
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
    imp_bfile_path='/oak/stanford/groups/mrivas/private_data/ukbb/24983/imp/pgen/ukb_imp_chr${SLURM_ARRAY_TASK_ID}_v2.mac1.hrc'
    cal_bfile_path='/oak/stanford/groups/mrivas/private_data/ukbb/24983/cal/pgen/ukb24983_cal_cALL_v2'
    # TODO: systematically name the output
    pheName=os.path.basename(pheFile).split('.')[0]
    outFile=os.path.join(outDir, 'ukb24983_v2.{0}.{1}'.format(pheName, kind))
    if kind == 'imputed':
        outFile += '.chr${SLURM_ARRAY_TASK_ID}'
        # make this an array job with nJobs=22
        cmd = make_plink_command(bpFile  = imp_bfile_path,
                                 pheFile = pheFile,
                                 outFile = outFile,
                                 pop     = pop,
                                 related = related,
                                 plink1  = plink1,
                                 arrayCovar = True)
    elif kind == 'genotyped': 
        # need one plink call with array a covariate, and one without
        outFile1 = outFile+'.both_arrays'
        cmd1 = make_plink_command(bpFile  = cal_bfile_path,
                                  pheFile = pheFile,
                                  outFile = outFile1,
                                  pop     = pop,
                                  related = related,
                                  plink1  = plink1,
                                  arrayCovar = True)
        outFile2 = outFile+'.one_array'
        cmd2 = make_plink_command(bpFile  = cal_bfile_path,
                                  pheFile = pheFile,
                                  outFile = outFile2,
                                  pop     = pop,
                                  related = related,
                                  plink1  = plink1,
                                  arrayCovar = False)
        # TODO: make sure that the plink output is actually getting formatted like this
        # join the plink calls, add some bash at the bottom to combine the output
        cmd = "\n\n".join([cmd1, cmd2] +  # this is the plink part, below joins the two files
                          ["if [ -f {0}.*.{3} ] cat {0}.*.{3} {1}.*.{3} | sort -k1,1n -k2,2n > {2}.{3}".format(
                               outFile1, outFile2, outFile, suffix) for suffix in ['glm.linear', 'glm.logistic.hybrid']] + 
                          ["cat {0}.log {1}.log > {2}.log".format(outFile1, outFile2, outFile),
                           "rm {0}.* {1}.*".format(outFile1, outFile2)])
    else:
        raise ValueError("argument kind must be one of (imputed, genotyped): {0} was provided".format(kind))
    sbatch = make_batch_file(batchFile = os.path.join(logDir, "gwas.{0}.{1}.sbatch.sh".format(kind,pheName)),
                             plinkCmd  = cmd,
                             memory    = memory,
                             time      = time,
                             partitions = partition)
    os.system(" ".join(("sbatch", "--array=1{}".format("-22" if kind == 'imputed' else ""), sbatch))) 


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
    parser.add_argument('--population', dest="pop", required=False, default=["white_british"], nargs=1,
                            help='Flag to indicate which ethnic group to use for GWAS. Must be one of all, white_british, e_asian, s_asian, african')
    parser.add_argument('--keep-related', dest="relatives", action='store_true',
                            help='Flag to keep related individuals in GWAS. Default is to remove them.')
    parser.add_argument('--batch-memory', dest="sb_mem", required=False, default=["16000"], nargs=1,
                            help='For underlying batch job submission: Amount of memory (in MB) to request. Default is 16000.')
    parser.add_argument('--batch-time', dest="sb_time", required=False, default=["24:00:00"], nargs=1,
                            help='For underlying batch job submission: Amount of time (DD-HH:MM:SS) to request. Default is 24 hours.')
    parser.add_argument('--batch-partitions', dest="sb_parti", required=False, default=["normal","owners"], nargs='*',
                            help='For underlying batch job submission: Compute partition to submit jobs. Default is normal,owners.')
    parser.add_argument('--log-dir', dest="log", required=False, nargs=1, default=[''],
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
    
 

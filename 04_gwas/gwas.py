#!/usr/bin/env python3
import os
from random import randint
_README_='''
A script for running GWAS with UK Biobank data (array/exome/imputed/cnv genotypes) using PLINK.

Author: Matthew Aguirre (SUNET: magu)

(Updated 10/14/2019 by E Flynn for sex-div analysis)
(Updated 11/06/2019 by Yosuke Tanigawa; add support for array_combined and array_imp_combined; code clean up.)
(Updated 11/27/2019 by Yosuke Tanigawa; add --keep option; code clean up.)
'''

def filterPopFile(outDir, popFile='', keepSexFile=''):
            # reformat the popFile to only contain the iids from one sex
        if not popFile or popFile=='':
            popFileSexFiltered=keepSexFile
        else:
            pop_iids = []
            keep_iids = []
            with open(popFile, 'r') as pop_f:
                for line in pop_f:
                    iid1, iid2 = line.split()
                    pop_iids.append(iid1.strip())
            with open(keepSexFile, 'r') as sex_f:
                for line in sex_f:
                    iid1 = line.strip()
                    keep_iids.append(iid1)
            filt_iids = set(pop_iids).intersection(set(keep_iids))
            # TODO: update to use temporary file (but that remains open thru the plink run?)
            popFileSexFiltered='{0}/tmp_pop_sex_{1}.keep'.format(outDir, randint(10000,99999))
            with open(popFileSexFiltered, 'w') as keep_f:
                for filt_iid in filt_iids:
                    keep_f.write("{0}\t{1}\n".format(filt_iid, filt_iid))
        return popFileSexFiltered

def updateKeepFile(outDir, qcDir, keepFile=None, pop=None, sexDiv=False, keepSexFile=''):
    if((keepFile is None) and (pop is not None) and (pop != 'all')):
        keepFile=os.path.join(qcDir, 'population_stratification','ukb24983_{}.phe'.format(pop))
    if sexDiv:
        # filter the population file to only contain IIDs of the specified sex        
        keepFile=filterPopFile(outDir, keepFile, keepSexFile)
    return(keepFile)

def make_plink_command(bpFile, pheFile, outFile, outDir, pop, keepFile=None, cores=None, memory=None, related=False, plink1=False, 
                       variantSubsetStr='', arrayCovar=False, sexDiv=False, keepSex='', keepSexFile='', includeX=False, maf=None, rmadd='', plink_opts=''):
    # paths to plink genotypes, input phenotypes, output directory are passed
    qcDir         = '/oak/stanford/groups/mrivas/ukbb24983/sqc/'
    keepFile=updateKeepFile(outDir, qcDir, keepFile=keepFile, pop=pop, sexDiv=sexDiv, keepSexFile=keepSexFile)
    unrelatedFile = os.path.join(qcDir,'ukb24983_v2.not_used_in_pca.phe') if not related else ''
    is_cnv_burden = (os.path.basename(bpFile[1]) == 'burden') or (os.path.basename(bpFile[1]) == 'ukb24983_cal_hla_cnv') or (os.path.basename(bpFile[1]) == 'ukb24983_hg19_cal_hla_cnv_imp')
    if os.path.dirname(pheFile).split('/')[-1] == "binary": 
        is_biomarker_binary = True
    elif len(os.path.basename(pheFile).split('_')) < 2:
        is_biomarker_binary = False
    else:
        is_biomarker_binary = os.path.basename(pheFile).split('_')[1]  == 'binary'
    if is_biomarker_binary:
        covarFile = '/oak/stanford/groups/mrivas/projects/biomarkers/covariate_corrected/output/covariates/logistic.covariates.phe'
    else:
        covarFile = os.path.join(qcDir, 'ukb24983_GWAS_covar.phe')    

    if (plink1 or (bpFile[0] == 'bfile')):
        genotypeStr='--bfile {}'.format(bpFile[1])
    elif (bpFile[0] == 'bpfile'):
        genotypeStr='--bpfile {}'.format(bpFile[1])
    elif (bpFile[0] == 'pfile_vzs'):
        genotypeStr='--pfile {} vzs'.format(bpFile[1])
    elif (bpFile[0] == 'pfile'):
        genotypeStr='--pfile {}'.format(bpFile[1])
    else:
        raise ValueError("Error: unsupported genotype file flag ({0})".format(bpFile[0]))
        
    # paste together the command from constituent parts
    if unrelatedFile and len(rmadd) > 0:
        os.system('cat ' + unrelatedFile + ' ' + rmadd + ' > tmp')
    cmd_plink = " ".join([
        "plink" if plink1 else "plink2",
        "--threads {0}".format(cores) if (cores is not None) else "",
        "--memory {0}".format(memory) if (memory is not None) else "",
        genotypeStr,
        "--chr 1-22" + (",X,XY" if includeX else ""),
        "--maf {0}".format(maf) if (maf is not None) else "",
        "--pheno", pheFile, "--pheno-quantile-normalize",
        "--glm firth-fallback hide-covar omit-ref ", "no-x-sex" if includeX else "",
        "--keep {0}".format(keepFile) if (keepFile is not None) else '', 
        "--remove {0}".format(unrelatedFile) if unrelatedFile and len(rmadd) == 0 else "",
        "--remove {0}".format(rmadd) if len(rmadd) > 0 and not unrelatedFile else "",
        "--remove tmp" if len(rmadd) > 0 and unrelatedFile else "",
        variantSubsetStr,
        "--covar", covarFile, 
        "--covar-name age ", "sex " if not sexDiv else "", 
        "Array " if arrayCovar else "", 
        "PC1-PC10", 
        "N_CNV LEN_CNV" if is_cnv_burden else "FastingTime" if is_biomarker_binary else "",
        "--covar-variance-standardize",
        "--vif 100000000" if pop in ['non_british_white', 'african', 'e_asian', 's_asian'] else "",
        "--out", outFile,
        plink_opts
    ])
    # gwas_sh=os.path.join(os.path.dirname(__file__), '04_gwas_misc.sh')
    gwas_sh="/oak/stanford/groups/mrivas/users/guhan/repos/ukbb-tools/04_gwas/04_gwas_misc.sh"
#    print(os.path.dirname(__file__))
    cmds = [
        "source {0}".format(gwas_sh),
        cmd_plink,
        "post_processing {0}".format(outFile),
    ]
    return("\n\n".join(cmds))

def make_plink_commands_arrayCovar(bpFile, outFile, make_plink_command_common_args, cores):
        # needs one plink call with genotyping array a covariate, and one without
        one_array_variants_txt='/oak/stanford/groups/mrivas/private_data/ukbb/24983/sqc/one_array_variants.txt'
        outFile1 = outFile+'.both_arrays'
        cmd1 = make_plink_command(
            bpFile  = bpFile,
            outFile = outFile1,
            arrayCovar = True,
            variantSubsetStr = "--exclude {0}".format(one_array_variants_txt),
            **make_plink_command_common_args
        )
        outFile2 = outFile+'.one_array'
        cmd2 = make_plink_command(
            bpFile  = bpFile,
            outFile = outFile2,
            arrayCovar = False,
            variantSubsetStr = "--extract {0}".format(one_array_variants_txt),
            **make_plink_command_common_args
        )
        if(cores is None):
            cores=1
        # join the plink calls, add some bash at the bottom to combine the output  
        return("\n\n".join([
            cmd2, cmd1,
            "combine_two_sumstats {0} {1} {2} {3}".format(outFile1, outFile2, outFile, cores)
        ]))

def make_batch_file(batchFile, plinkCmd, cores, memory, time, partitions):
    with open(batchFile, 'w') as f:
       # formats options for sbatch header and pastes input command below it
       f.write("\n".join([
           "#!/bin/bash","",
           "#SBATCH --job-name=RL_GWAS",
           "#SBATCH --output={}".format(os.path.join(os.path.dirname(batchFile), "rl-gwas.%A-%a.out")),
           "#SBATCH  --error={}".format(os.path.join(os.path.dirname(batchFile), "rl-gwas.%A-%a.err")),
           "#SBATCH --cores={}".format(cores),
           "#SBATCH --mem={}".format(memory),
           "#SBATCH --time={}".format(time),
           "#SBATCH -p {}".format(','.join(partitions)),
           '#SBATCH --constraint="CPU_GEN:HSW|CPU_GEN:BDW|CPU_GEN:SKX"', # plink2 avx2 compatibility
           'set -beEuo pipefail',
           "", plinkCmd
       ]))
    return(batchFile)


def run_gwas(kind, pheFile, outDir='', pop='white_british', keepFile=None, related=False, plink1=False, 
             logDir=None, cores="4", memory="24000", rmAdd=None, time="1-00:00:00", partition=["normal","owners"], now=False,
             sexDiv=False, keepSex='', keepSexFile='', includeX=False, plink_opts=''):
    # ensure usage
    rmadd = ""
    if not os.path.isfile(pheFile):
        raise ValueError("Error: phenotype file {0} does not exist!".format(pheFile))
    if rmAdd is not None:
        rmadd = rmAdd
    if not os.path.isdir(outDir):
        raise ValueError("output directory {0} does not exist!".format(outDir))
    if ((logDir is None) or (not os.path.isdir(logDir))):
        logDir=os.getcwd()
        print(f'Warning: Logging directory either does not exist or was unspecified. Defaulting to {logDir}...')
    if sexDiv and not os.path.isfile(keepSexFile):
        raise ValueError("Error: keep sex file {0} does not exist!".format(keepSexFile))
    if sexDiv and kind != "genotyped":
        print("Warning: sex div GWAS has only been tested on array data! Please make sure to test this more or use --run-array.")
    if pop not in ['all', 'white_british', 'non_british_white', 'african', 's_asian', 'e_asian']:
        raise ValueError("population must be one of (all, white_british, non_british_white, african, s_asian, e_asian)")
    # paths for running gwas
    pgen_root='/oak/stanford/groups/mrivas/private_data/ukbb/24983/'
    genotype_file={
        'imputed':            ('pfile_vzs', os.path.join(pgen_root,'imp','pgen','ukb24983_imp_chr${SLURM_ARRAY_TASK_ID}_v3')),
        'genotyped':          ('bpfile',    os.path.join(pgen_root,'cal','pgen','ukb24983_cal_cALL_v2_hg19')),
        'array-combined':     ('bpfile',    os.path.join(pgen_root,'array_combined','pgen','ukb24983_cal_hla_cnv')),
        'array-imp-combined': ('pfile_vzs', os.path.join(pgen_root,'array_imp_combined','pgen','ukb24983_hg19_cal_hla_cnv_imp')),
        'cnv':                ('bpfile',    os.path.join(pgen_root,'cnv','pgen','cnv') + ' --mac 15'),
        'cnv-burden':         ('bpfile',    os.path.join(pgen_root,'cnv','pgen','burden')),
        'exome-spb':          ('bpfile',    os.path.join(pgen_root,'exome','pgen','spb','data','ukb_exm_spb')),
        'exome-fe':           ('bpfile',    os.path.join(pgen_root,'exome','pgen','fe','data','ukb_exm_fe')),
        'hla':                ('bpfile',    os.path.join(pgen_root,'hla','pgen','ukb_hla_v3'))
    }
    pheName=os.path.basename(pheFile).split('.')[0]
    outFile=os.path.join(outDir, 'ukb24983_v2_{0}.{1}{2}.{3}'.format(
        'hg38' if 'exome' in kind else 'hg19', 
        pheName, 
        "_{0}".format(keepSex) if sexDiv else "", 
        kind
    ))
    
    make_plink_command_common_args = {
        'pheFile' : pheFile,
        'outDir'  : outDir,
        'pop'     : pop,
        'related' : related,
        'plink1'  : plink1,
        'cores'   : cores,
        'memory'  : memory,
        'sexDiv'  : sexDiv,
        'keepFile' : keepFile,
        'keepSex' : keepSex,
        'keepSexFile' : keepSexFile,
        'includeX' : includeX, 
        'rmadd' : rmadd,
        'plink_opts': plink_opts
    }
    
    # this is where the fun happens
    if kind == 'imputed':
        # needs one array job per chromosome (for compute time), hence the slurm variable
        outFile += '.chr${SLURM_ARRAY_TASK_ID}'
        cmd = make_plink_command(
            bpFile  = genotype_file[kind],
            outFile = outFile,
            maf = 0.005,
            arrayCovar = True, # why we have array for the imputation?
            **make_plink_command_common_args
        )
    elif kind in ['genotyped', 'array-combined', 'array-imp-combined']:
        cmd = make_plink_commands_arrayCovar(
            genotype_file[kind], outFile, make_plink_command_common_args, cores
        )
    # more usage management, in case someone wants to import the function for use elsewhere
    elif kind in ['cnv', 'cnv-burden', 'exome-spb', 'exome-fe', 'hla']:
        cmd = make_plink_command(
            bpFile  = genotype_file[kind],
            outFile = outFile,
            **make_plink_command_common_args
        ) 
    else:
        raise ValueError("argument kind must be one of (imputed, genotyped, exome-spb, exome-fe): {0} was provided".format(kind))
    # run immediately, OR make the batch job submission file, then call it with an appropriate array handler
    if now:
        now_cmd = '\n\n'.join(['set -beEuo pipefail', cmd])
        print("Running the below: \n'''\n" + now_cmd + "\n'''\n") 
        os.system(now_cmd)
    else:
        sbatch = make_batch_file(batchFile = os.path.join(logDir, "gwas.{0}.{1}.sbatch.sh".format(kind,pheName)),
                                 plinkCmd  = cmd,
                                 cores     = cores,
                                 memory    = memory,
                                 time      = time,
                                 partitions = partition)
        os.system(" ".join(("sbatch", "--array=1{}".format("-22" if kind == 'imputed' else ""), sbatch))) 
    return


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
    parser.add_argument('--run-exome', dest="ex1", action='store_true',
                            help="Run GWAS on exome data (Regeneron calls)")
    parser.add_argument('--run-exome-gatk', dest="ex2", action='store_true',
                            help="Run GWAS on exome data (GATK calls)")
    parser.add_argument('--run-imputed', dest="imp", action='store_true',
                            help='Run GWAS on imputed data') 
    parser.add_argument('--run-cnv', dest="cnva", action='store_true',
                            help='Run GWAS on array-derived CNV genotypes') 
    parser.add_argument('--run-cnv-burden', dest="cnvb", action='store_true',
                            help='Run CNV burden test (GWAS on 0/1 CNV overlaps gene, from array-derived CNV genotypes)') 
    parser.add_argument('--run-hla', dest="hla", action='store_true',
                            help='Run GWAS on imputed HLA allelotypes') 
    parser.add_argument('--run-array-imp-combined', dest="array_imp_combined", action='store_true',
                            help='Run GWAS on the array_imp_combined dataset')
    parser.add_argument('--run-array-combined', dest="array_combined", action='store_true',
                            help='Run GWAS on the array_combined dataset')
    parser.add_argument('--pheno', dest="pheno", required=True,
                            help='Path to a phenotype file')
    parser.add_argument('--out', dest="outDir", required=True, 
                            help='Path to desired output *directory*. Summary stats will be output according to phenotype name (derived from passed file) and Rivas Lab specification for GBE. Defaults to current working directory.')
    parser.add_argument('--population', dest="pop", required=False, default="white_british",
                            help='Flag to indicate which ethnic group to use for GWAS. Must be one of all, white_british, non_british_white, e_asian, s_asian, african')
    parser.add_argument('--keep', dest="keep", required=False, default=None,
                            help='A file that specifies the list of individuals for GWAS. It can be overwritten by sex-specific subcommands')
    parser.add_argument('--remove-add', dest="rmadd", required=False, default=None,
                            help='Flag to indicate file of list of additional individuals to remove')
    parser.add_argument('--keep-related', dest="relatives", action='store_true',
                            help='Flag to keep related individuals in GWAS. Default is to remove them.')
    parser.add_argument('--cores', dest="cores", required=False, default=4, type=int,
                            help='For underlying plink command/batch job submission: Amount of cores to request. Default is 4 cores for batch jobs.')    
    parser.add_argument('--memory', dest="mem", required=False, default=24000, type=int,
                            help='For underlying plink command/batch job submission: Amount of memory (in MB) to request. Default is 24000 for batch jobs.')
    parser.add_argument('--batch-time', dest="sb_time", required=False, default="24:00:00",
                            help='For underlying batch job submission: Amount of time (DD-HH:MM:SS) to request. Default is 24 hours.')
    parser.add_argument('--batch-partitions', dest="sb_parti", required=False, default=["normal","owners"], nargs='*',
                            help='For underlying batch job submission: Compute partition to submit jobs. Default is normal,owners.')
    parser.add_argument('--log-dir', dest="log", required=False, default=None,
                            help="Directory in which to place log files from this script and its related SLURM jobs. Default is current working directory.")
    parser.add_argument('--run-now', dest="local", action='store_true',
                            help='Flag to run GWAS immediately, on the local machine, rather than submitting a script with sbatch. Not available for use with --run-imputed.')
    parser.add_argument('--sex-div', dest="sex_div",required=False, default=False,
                            help='Whether to run sex-divided GWAS, this removes sex as a covariate from the analysis. Default is False. If specifying this, please also use --keep-sex to specify the sex to keep and --keep-sex-file to specify the file location.')
    parser.add_argument('--keep-sex', dest="keep_sex", required=False, default='',
                            help='Which sex to keep, will be used in constructing the output filenames, for use with --sex-div.')
    parser.add_argument('--keep-sex-file', dest="keep_sex_file", required=False, default='',
                            help='Location of the file specifying the IIDs to include related to that sex, for use with --sex-div.')
    parser.add_argument('--include-x', dest="include_x", action='store_true', default=False,
                            help='Whether to include the X chromosome, defaults to False')
    parser.add_argument('--additional-plink-opts', dest="plink_opts", required=False, default=[], nargs='*',
                            help='Addtional options for plink')

    args = parser.parse_args()
    # TODO: feature add: genotype model  
    print(args) 
    # ensure handler-relevant usage (more insurance is in run_gwas()):
    flags = [args.imp, args.arr, args.ex1, args.ex2, args.cnva, args.cnvb, args.hla, args.array_imp_combined, args.array_combined]
    kinds = ['imputed','genotyped','exome-spb','exome-fe','cnv','cnv-burden', 'hla', 'array-imp-combined', 'array-combined']
    if not any(flags):
        raise ValueError("Error: no analysis specified, did you mean to add --run-array?")
    if args.local and args.imp:
        raise ValueError("--run-imputed cannot be present in conjunction with --run-now!")
    if (args.sex_div and args.keep_sex == '') or (args.sex_div and args.keep_sex_file == ''):
        raise ValueError("Sex-div analysis is indicated but either the sex to keep or the file for that is missing. Please use --keep-sex and --keep-sex-file to specify both.")
    for flag,kind in filter(lambda x:x[0], zip(flags,kinds)):
        run_gwas(kind=kind, pheFile=os.path.realpath(args.pheno), outDir=args.outDir,                 
                 pop=args.pop, rmAdd = args.rmadd, keepFile=args.keep, related=args.relatives, plink1=args.plink1, 
                 logDir=args.log, cores=args.cores, memory=args.mem,
                 time=args.sb_time, partition=args.sb_parti, now=args.local, 
                 sexDiv=args.sex_div, keepSex=args.keep_sex, keepSexFile=args.keep_sex_file,
                 includeX=args.include_x, plink_opts=' '.join([f'--{x}' for x in args.plink_opts]))

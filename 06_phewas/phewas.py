#!/bin/bash
import os
import glob

# useful
os.system('ml load plink2')
_README_ = '''
A script to run phewas. More documentation forthcoming.

Author: Matthew Aguirre (SUNET: magu)
'''

# Input: gene name(s) / variant names(s) / region
# variant names for imputed data are chrom:pos:ref:alt

cal_bims = glob.glob('/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_chr*_v2.bim')
imp_bims = glob.glob('/oak/stanford/groups/mrivas/ukbb24983/imp/pgen/ukb_imp_chr*_v2.mac1.hrc.bim')
wex_bims = ['/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/ukb24983_exome.bim']

# 1. Get all variants in gene (if applicable)
def find_variants(chrom,start,end,exome=False):
    var_list = {} #file:[list of variant ids]
    # helper
    overlaps = lambda x,c,bp1,bp2: str(c)==x[0] and int(bp1) <= int(x[3]) and int(x[3]) <= int(bp2)
    bims = wex_bims if exome else filter(lambda x:'_chr{}_'.format(chrom) in x, cal_bims+imp_bims)
    for bim in bims:
        with open(bim, 'r') as f:
            var_list[bim[:-4]] = [line.split()[1] for line in f if overlaps(line.split(),chrom,start,end)] 
    return var_list

def find_named_variants(var_list,exome=False):
    file_map = {} # this is the same as var_list in find_variants()
    for bim in wex_bims if exome else cal_bims + imp_bims:
        with open(bim, 'r') as f:
            # find variants, update list of variants to find
            file_map[bim[:-4]] = [line.split()[1] for line in f if line.split()[1] in var_list]
            var_list = list(set(var_list).difference(file_map[bim[:-4]])) 
    if len(var_list) > 0:
        print("Could not find the following variants:\n{}".format("\n".join(var_list)))
        print("Note: variant names for imputed/exome data are formatted CHR:POS_REF_ALT, and exomes are in hg38!")
    return {k:v for k,v in file_map.items() if len(v) > 0}    

            
def get_variants(gene,exome=False):
    with open('genes_hg{}.txt'.format(38 if exome else 19), 'r') as f:
        for line in f:
            if gene == line.rstrip().split()[-1]:
                c,bp1,bp2 = line.split()[:3]
                return find_variants(c,bp1,bp2,exome)
    raise ValueError("Could not find gene: {}".format(gene))

# 2. Run PheWAS(es)
def run_phewas(bfile, var_ids, indfile=None, norm=True, out_prefix='phewas'):
    master_phe = '/oak/stanford/groups/mrivas/dev-ukbb-tools/phewas/resources/master.phe'
    covars     = '/oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar.phe'
    # write out temp file for variants
    with open(out_prefix + '.varlist.txt', 'w') as o:
        o.write('\n'.join(var_ids))
    # do the thing
    os.system(' '.join(['plink2 --bfile', bfile, '--pheno', master_phe, 
                               '--pheno-quantile-normalize --covar-variance-standardize ' if norm else '',
                               '--out', out_prefix, 
                               '--keep {}'.format(indfile) if indfile is not None else '',
                               '--covar', covars, '--covar-name age sex PC1-PC4',
                               '--extract', out_prefix + '.varlist.txt', 
                               '--glm firth-fallback hide-covar']))
    return
     

# 3. Combine results
def combine_output(common_prefix, out_prefix, keep_na=False):
    # do qts and bins separately
    for is_bin, ending in enumerate(['linear', 'logistic.hybrid']):
        files = glob.glob(common_prefix + '*' + ending)
        print(len(files))
        with open(out_prefix + 'phewas.glm.' + ending, 'w') as o:
            for n,f in enumerate(files):
                pheno = os.path.basename(f).split('.')[-3 - is_bin]
                with open(f, 'r') as i:
                    for j,line in enumerate(i):
                        if n == 0 and j == 0:
                            o.write('#PHENO\t' + line[1:])
                        elif j > 0 and (keep_na or line.split()[-1] != 'NA'):
                            o.write(pheno + '\t' + line)
    return


# ok now let's actually do the darn thing:
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser( formatter_class=argparse.RawDescriptionHelpFormatter,
                                          description=_README_ )
    parser.add_argument('--gene', dest="gene", required=False, default = [], nargs=1,
                            help='input gene for phewas')
    parser.add_argument('--region', dest="cpra", required=False, default = [], nargs=1,
                            help='input region for phewas (format: CHROM:BP1-BP2)')
    parser.add_argument('--keep', dest="inds", required=False, default = [], nargs=1,
                            help='path to file containing individuals to keep (this takes precedence over --population)')
    parser.add_argument('--population', dest="pop", required=False, default = ['all'], nargs=1,
                            help='ethnic group to subset to for analysis: must be one of white_british,african,e_asian,s_asian,all')
    parser.add_argument('--variants', dest="vars", required=False, default = [], nargs='*', 
                            help='input variant ID(s) for phewas (can also be a file, with one variant ID per line)')
    parser.add_argument('--exome', dest="wex", action="store_true", default = False,  
                            help='flag to run phewas using exome data (works with all variant options, but note this data uses hg38!)')
    parser.add_argument('--qt-norm', dest="norm", action="store_true", default = False,  
                             help='flag to normalize quantitative phenotypes (equivalent to plink --pheno-quantile-normalize)')
    parser.add_argument('--out', dest="out", required=True, default = ['phewas'], nargs=1,
                             help='path to output (prefix, will be passed directly to plink)')
    args = parser.parse_args()
    print(args)
    
    # ensure usage
    if sum(map(lambda x: not x, (args.gene, args.cpra, args.vars))) < 2:
        raise ValueError("Only one of --gene, --region, --variants, can be supplied.")
    # get individuals
    if args.pop[0] == 'all':
        indfile = None
    elif args.inds:
        indfile = args.inds[0]
    elif args.pop[0] in ['white_british','african','e_asian','s_asian']:
        indfile = '/oak/stanford/groups/mrivas/private_data/ukbb/24983/sqc/population_stratification/ukb24983_'+args.pop[0]+'.phe'
    else:
        raise ValueError("--population must be one of white_british,african,e_asian,s_asian,all")
    # get variants
    if args.gene:
        bfile_to_vars = get_variants(args.gene[0],args.wex)
    elif args.cpra:
        chrom, bps = args.cpra[0].split(':')
        bp1, bp2   = bps.split('-')
        bfile_to_vars = find_variants(chrom,bp1,bp2,args.wex)
    elif args.vars:
        if os.path.isfile(args.vars[0]):
            with open(args.vars[0], 'r') as f:
                vs = [line.rstrip().split()[0] for line in f]
            print("{} variant IDs found in input file.".format(len(vs)))
        else:
            vs = args.vars
        bfile_to_vars = find_named_variants(vs,args.wex)
    else:
        raise ValueError("One of --gene, --region, --variants, must be supplied!")
    # run phewas
    print("{} variants found from specified input.".format(sum(map(len, bfile_to_vars.values()))))

    if sum(map(len,bfile_to_vars.values())) >= 100:
        print("WARNING: attempting phewas with > 100 variants! This will probably take awhile, and might not finish...")

    print("Running phewas...")
    
    # give these runs a random hash to prevent wrong files from being included when we join results
    import random
    temp_out = args.out[0] + '.' + str(int(1000000 * random.random()))
    for n,(bfile,variants) in enumerate(filter(lambda x: len(x[1]) > 0, bfile_to_vars.items())):
        run_phewas(bfile, variants, indfile, args.norm, temp_out+'.'+str(n))
    
    # now combine those files into the final result
    combine_output(temp_out, args.out[0], keep_na = False)
    
    # and remove the middling files
    os.system("rm {}*".format(temp_out))
    # done    

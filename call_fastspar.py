import os
import json
import argparse
from typing import List, Set, Dict, Tuple

import pandas as pd
from myutils import *

#Based on https://github.com/scwatts/fastspar

SEED = 123

def read_metadata(fname):
    if os.path.isfile(fname):
        metad = pd.read_csv(fname, sep='\t')    
    else:
        logger.log("WARNING: Metadata file not found!", bcolors.WARNING)
        metad = pd.DataFrame()
    return metad

def read_otu_table(fname):
    abndf = pd.read_csv(fname, sep='\t')
    abndf = abndf.set_index(abndf.columns[0], drop=True)
    abndf = abndf.select_dtypes(include='number')
    return abndf

def split_data(abundances, metadata, splitvars):
    if metadata.shape[0] == 0 or splitvars == '':
        logger.log("WARNING: Not splitting data into groups", bcolors.WARNING)
        return {"all":abundances}
    assert all([i in abundances.columns for i in metadata.sampleID.tolist()]), "Not all samples in Abundances table"
    splitvars = splitvars.split(',')
    splitmeta = metadata.groupby(splitvars)
    logger.log(f"WARNING: Splitting data into {len(splitmeta)} groups", bcolors.WARNING)
    res = {}
    for pname, pdf in splitmeta:
        pname = '_'.join(pname)
        filtab = abundances[pdf.sampleID]
        res.setdefault(pname, filtab)
    return res


def write_split_data(splitdata, original_name, outdir):
    fnames = []
    for name, df in splitdata.items():
        fname = f"{outdir}/splitAbn_{original_name}_{name}.csv"
        df=df.reset_index()
        df=df.rename(columns={df.columns[0]: "#OTU ID"})
        df.to_csv(fname, sep="\t", index = False)
        logger.log(f"Written split data: {fname}", bcolors.OKBLUE)
        fnames.append(fname)
    return fnames

def get_fastspar_commands(fname, outdir, fastspar_args, cleanup = True):
    sname = os.path.basename(fname)
    sname = os.path.splitext(sname)[0]
    btdir = f"{outdir}/bootstrap_counts_{sname}"
    btcordir = f"{outdir}/bootstrap_correlation_{sname}"

    correlation=f"{outdir}/median_correlation_{sname}.tsv"
    covariance=f"{outdir}/median_covariance_{sname}.tsv"
    pvalues=f"{outdir}/pvalues_{sname}.tsv"

    create_directories([btdir, btcordir])
    
    cmd1 = (f"fastspar --yes --otu_table {fname} "
            f"--correlation {correlation} "
            f"--covariance {covariance} "
            f"--iterations {fastspar_args['iterations']} "
            f"-x {fastspar_args['exclusion_iterations']} "
            f"-e {fastspar_args['exclusion_threshold']} "
            f"-t {fastspar_args['threads']} "           
            )
    
    cmd2 = (f"fastspar_bootstrap --otu_table {fname} "
            f"--number {fastspar_args['num_random']} "
            f"--prefix {btdir}/{sname} "
            f"-s {fastspar_args['seed']} ")
    
    cmd3 = ("parallel fastspar --yes --otu_table {} "
            f"--correlation {btcordir}"
            "/cor_{/} "
            f"--covariance {btcordir}"
            "/cov_{/} "
            f"-i {fastspar_args['iterations_parallel']} ::: "
            f"{btdir}/*")
    
    cmd4 = (f"fastspar_pvalues --otu_table {fname} "
            f"--correlation {correlation} "
            f"--prefix {btcordir}/cor_{sname} "
            f"--permutations {fastspar_args['num_permuts']} "
            f"-t {fastspar_args['threads']} " 
            f"--outfile {pvalues}")
    
    cmdlist = [cmd1, cmd2, cmd3, cmd4]

    if cleanup:
        cmd5 = f"rm -rf {btdir} {btcordir}"
        cmdlist.append(cmd5)
    return cmdlist

    #



def run_fastspar_all(fnames, outdir, fastspar_args, cleanup = True):
    for f in fnames:
        logger.log(f"Running FastSpar for file: {f}", bcolors.BOLD)
        cmds = get_fastspar_commands(f, outdir, fastspar_args, cleanup)
        for c in cmds:
            run_command(c)
        logger.log(f"Finished FastSpar for file: {f}", bcolors.OKBLUE)
    


# ----- command line parsing -----
parser: argparse.ArgumentParser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, 
                                                          description='Call FastSpar Wrapper')

group1 = parser.add_argument_group('General input and output options')

group1.add_argument('-m', '--metadata', type = str, help='Metadata with sampleID column', default="")
group1.add_argument('-a', '--abundances', type = str, help='Raw abundances (otu table)', default="")
group1.add_argument('-s', '--splitvars', type = str, help='Comma separated variables to split data into', default="")
group1.add_argument('-o', '--outdir', type = str, help='Output directory', default="FastSpar")
group1.add_argument('-c', '--cleanup', type = str, 
                    choices=['True', 'False', 'T', 'F'], help='Clean intermediate files', default='True')

group2 = parser.add_argument_group('FastSpar options')
group2.add_argument('-w', '--seed', type = int, help='Random seed', default=SEED)
group2.add_argument('-i', '--iterations', type = int, help='Number of iterations FastSpar', default=50)
group2.add_argument('-x', '--exclusion_iterations', type = int, help='Number of exclusion iterations FastSpar', default=10)
group2.add_argument('-j', '--iterations_parallel', type = int, help='Number of iterations from bootstrapped samples', default=5)
group2.add_argument('-n', '--nrand', type = int, help='Random samples', default=1000)
group2.add_argument('-p', '--permutations', type = int, help='Number of permutations', default=1000)
group2.add_argument('-e', '--exclusion_threshold', type = float, help='Correlation strength exclusion threshold.', default=0.1)
group2.add_argument('-t', '--threads', type = int, help='Threads for FastSpar first call and fastspar_pvalues.', default=1)  

 
# Example:
# python call_fastspar.py -o test -m /home/carmoma/Documents/EarlyCause/results/descriptive1/METADATA_FILT.tsv -a /home/carmoma/Documents/EarlyCause/results/descriptive1/OTUs_FILT.tsv -s Time,Group
def main():
    args = parser.parse_args()
    metadatafile: str = args.metadata
    abnfile: str = args.abundances
    splitvars: str = args.splitvars
    outdir: str = args.outdir
    cleanup: str = True if args.cleanup in ['T', 'True'] else False

    if cleanup:
        logger.log("WARNING: Cleaning intermediate files!", bcolors.WARNING)

    try:
        os.mkdir(outdir)
        logger.log(f"Directory created: {outdir}", bcolors.OKBLUE)
    except:
        logger.log(f"Directory creation error: {outdir}", bcolors.FAIL)
    
    fastspar_args = {}
    fastspar_args["seed"] = args.seed
    fastspar_args["num_random"] = args.nrand
    fastspar_args["num_permuts"] = args.permutations
    fastspar_args["iterations"] = args.iterations
    fastspar_args["exclusion_iterations"] = args.exclusion_iterations
    fastspar_args["iterations_parallel"] = args.iterations_parallel
    fastspar_args["exclusion_threshold"] = args.exclusion_threshold
    fastspar_args["threads"] = args.threads

    with open(f"{outdir}/fastspar_options.json", "w") as outfile: 
        json.dump(fastspar_args, outfile)
        logger.log(f"Printed FastSpar options: {outdir}/fastspar_options.json", bcolors.OKBLUE)

    metadata: pd.DataFrame = read_metadata(metadatafile)
    abundances: pd.DataFrame = read_otu_table(abnfile)

    splitdf: List[pd.DataFrame] = split_data(abundances, metadata, splitvars)

    original_name: str = os.path.basename(abnfile).replace('.tsv', '')

    fnames: List[str] = write_split_data(splitdf, original_name, outdir)

    run_fastspar_all(fnames, outdir, fastspar_args, cleanup)


    


if __name__ == "__main__":
    main()
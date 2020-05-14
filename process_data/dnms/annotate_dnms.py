import os
import yaml
import glob
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-y", "--config",
                    help="config file",
                    nargs='?',
                    metavar='',
                    type=str)

parser.add_argument("-t", "--test",
                    help="use test vcf",
                    action="store_true")

args = parser.parse_args()

with open(args.config, "r") as f:
    config = yaml.load(f)


###############################################################################
# print to stderr
###############################################################################
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

sourcedir = config['sourcedir']
outdir = config['outdir']

###############################################################################
# create output path if it doesn't already exist
###############################################################################
# outdir = os.path.realpath(args.outdir)

if not os.path.exists(outdir):
    os.makedirs(outdir)

if args.test:
    vcfdir = sourcedir + "sandbox_data/"
    vcftail = ".freeze3a.test.vcf"
else:
    vcfdir = config['vcfdir']
    vcftail = ".freeze5.dnms.vcf"

pythonpath = config['pythonpath']
bcftoolspath = config['bcftoolspath']

vcfs = glob.glob(vcfdir)

for i in range(1,23):
    chrom = "chr" + str(i)
    m2r = "bcftools +missing2ref " + vcfdir + chrom + vcftail
    annocmd = pythonpath + " " + sourcedir + "scripts/annotate_motifs.py -s -i - -f " + config['grch38'] + " " + "-l 7"
    
    querycols = "'%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/type\\t%INFO/motif\\t%INFO/sample\\n'"
    querycmd = bcftoolspath + " query -f " + querycols
    sortcmd = "sort -k2,2n"
    outcmd = " > " + outdir + chrom + ".freeze5.dnms.txt"
    
    cmd = m2r + " | " + annocmd + " | " + querycmd + " | " + sortcmd + outcmd
    eprint(cmd)
    os.system(cmd)

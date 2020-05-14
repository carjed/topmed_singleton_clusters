import os
import yaml
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--chr",
                    help="chromosome to run",
                    nargs='?',
                    type=int,
                    metavar='',
                    default=1)

parser.add_argument("-b", "--begin",
                    help="start of region",
                    nargs='?',
                    type=int,
                    metavar='',
                    default=1)
                    
parser.add_argument("-e", "--end",
                    help="end of region",
                    nargs='?',
                    type=int,
                    metavar='',
                    default=1000000)

parser.add_argument("-s", "--samplefile",
                    help="file with sample IDs to include (one per line)",
                    nargs='?',
                    metavar='',
                    type=str)

parser.add_argument("-o", "--outdir",
                    help="directory to store output files \
                        (do NOT include a trailing '/')",
                    nargs='?',
                    type=str,
                    metavar='',
                    default="~/")

parser.add_argument("-m", "--mask",
                    help="filter sites passing mask",
                    action="store_true")

parser.add_argument("-t", "--test",
                    help="use test vcf",
                    action="store_true")

parser.add_argument("-C", "--config",
                    help="path to yaml config with necessary paths",
                    nargs='?',
                    type=str,
                    metavar='',
                    default="./config.yaml")

args = parser.parse_args()

with open(args.config, "r") as f:
    config = yaml.safe_load(f)

###############################################################################
# create output path if it doesn't already exist
###############################################################################
outdir = os.path.realpath(args.outdir)

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

###############################################################################
# print to stderr
###############################################################################
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

###############################################################################
# process args
###############################################################################
chrom = "chr" + str(args.chr)

bcfregion = chrom + ":" + str(args.begin) + "-" + str(args.end)
outregion = bcfregion.replace(":", ".")

filtercmd = config["bcftoolspath"] + " view -f PASS -Ou -S " + args.samplefile + " -r " + bcfregion + " " + config["vcfdir"] + config["vcfhead"] + chrom + config["vcftail"]
viewcmd = config["bcftoolspath"] + " view -c 1 -C 1 -Ou" 
maskcmd = config["bcftoolspath"] + " view -T " + config["maskfile"] + " -Ou"
annocmd = config["pythonpath"] + " " + config["sourcedir"] + "process_data/annotate_motifs.py -s -i - -f " + config["grch38"] + " -l 7"

querycols = "'%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/type\\t%INFO/motif\\t%INFO/sample\\n'"
querycmd = config["bcftoolspath"] + " query -f " + querycols + " > " + args.outdir + "/" + outregion + ".freeze5.singletons.txt"

# headercmd = "/bin/sed '1s/.*/CHR\\tPOS\\tREF\\tALT\\tTYPE\\tMOTIF\\tID\\n&/'"
# distcmd = pythonpath + sourcedir + "scripts/add_dist.py -i - > /net/snowwhite/home/jedidiah/" + chrom + ".freeze3.singletons.txt"

# less -NS" #
# cmd = "bcftools view -c 1 -C 1 " + vcfdir + chrom + vcftail + " | python annotate_motifs.py -i - -f ~/GRCh37-lite.fa -l 7 | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/type\t%INFO/motif\t%INFO/sample\n' | sed '1s/.*/CHR\tPOS\tREF\tALT\tTYPE\tMOTIF\tID\n&/'  | python add_dist.py -i - > ~/chr1.freeze3.singletons.txt"

# run cmd
# cmd = filtercmd + " | " + viewcmd + " | " + annocmd + " | less -NS"
if args.mask:
    # viewcmd = bcftoolspath + "view -c 1 -C 1 -Ou" 
    cmd = filtercmd + " | " + viewcmd + " | " + maskcmd + " | " + annocmd + " | " + querycmd # + " | " + headercmd + " | " + distcmd
else:
    # viewcmd = bcftoolspath + "view -c 1 -C 1 -Ou" 
    cmd = filtercmd + " | " + viewcmd + " | " + annocmd + " | " + querycmd # + " | " + headercmd + " | " + distcmd
    
eprint(cmd)
os.system(cmd)

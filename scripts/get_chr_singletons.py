import os
# import yaml
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

parser.add_argument("-t", "--test",
                    help="use test vcf",
                    action="store_true")

args = parser.parse_args()

# with open(args.config, "r") as f:
#     config = yaml.load(f)

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

chrom = "chr" + str(args.chr)

sourcedir = "/net/snowwhite/home/jedidiah/"

if args.test:
    vcfdir = sourcedir + "sandbox_data/"
    vcftail = ".freeze3a.test.vcf"
else:
    vcfdir = "/net/topmed/working/dtaliun/TOPMed_paper_phased/unrelated/"
    vcftail = ".freeze3a.gtonly.Eagle.Phased.unrelated_samples_for_analysis.min_ac_1.vcf.gz"

pythonpath = sourcedir + "anaconda3/envs/anno/bin/python"
bcftoolspath = sourcedir + "anaconda3/envs/anno/bin/bcftools"

# vcftest = "~/chr9.freeze3a.test.vcf"

bcfregion = str(args.chr) + ":" + str(args.begin) + "-" + str(args.end)
# bcfregion = chrom + ":" + str(args.begin) + "-" + str(args.end)
outregion = bcfregion.replace(":", ".")

filtercmd = bcftoolspath + "view -Ou -S " + args.samplefile + " -r " + bcfregion + " " + vcfdir + chrom + vcftail
viewcmd = bcftoolspath + "view -c 1 -C 1 -Ov" 
annocmd = pythonpath + " " + sourcedir + "scripts/annotate_motifs.py -s -i - -f " + sourcedir + "GRCh37-lite.fa -l 7"

querycols = "'%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/type\\t%INFO/motif\\t%INFO/sample\\n'"
querycmd = bcftoolspath + "query -f " + querycols + " > " + args.outdir + "/chr" + outregion + ".freeze3.singletons.txt"

# headercmd = "/bin/sed '1s/.*/CHR\\tPOS\\tREF\\tALT\\tTYPE\\tMOTIF\\tID\\n&/'"
# distcmd = pythonpath + sourcedir + "scripts/add_dist.py -i - > /net/snowwhite/home/jedidiah/" + chrom + ".freeze3.singletons.txt"

# less -NS" #
# cmd = "bcftools view -c 1 -C 1 " + vcfdir + chrom + vcftail + " | python annotate_motifs.py -i - -f ~/GRCh37-lite.fa -l 7 | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/type\t%INFO/motif\t%INFO/sample\n' | sed '1s/.*/CHR\tPOS\tREF\tALT\tTYPE\tMOTIF\tID\n&/'  | python add_dist.py -i - > ~/chr1.freeze3.singletons.txt"

# run cmd
# cmd = filtercmd + " | " + viewcmd + " | " + annocmd + " | less -NS"
cmd = filtercmd + " | " + viewcmd + " | " + annocmd + " | " + querycmd # + " | " + headercmd + " | " + distcmd
eprint(cmd)
os.system(cmd)

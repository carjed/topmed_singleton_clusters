import os
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--chr",
                    help="chromosome to run",
                    nargs='?',
                    type=int,
                    metavar='',
                    default=1)

parser.add_argument("-p", "--pop",
                    help="population data to aggregate (either eur or afr)",
                    nargs='?',
                    metavar='',
                    type=str)

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
# print to stderr
###############################################################################
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


#-----------------------------------------------------------------------------
# define paths
#-----------------------------------------------------------------------------
chrom = "chr" + str(args.chr)

#-----------------------------------------------------------------------------
# concatenate 1Mb chunks from each chromosome
# uses 'ls -v' to list 1Mb chunk files in correct (numeric) order
#-----------------------------------------------------------------------------
# sourcedir + "topmed_freeze3_singletons/" + args.pop + "/" + \
concatcmd = "ls -v " + \
    config["sourcedir"] + "topmed_freeze5_singletons_mask/" + args.pop + "/" + \
    chrom + ".* | xargs cat"

#-----------------------------------------------------------------------------
# sed command to prepend header
#-----------------------------------------------------------------------------
headercmd = "/bin/sed '1s/.*/CHR\\tPOS\\tREF\\tALT\\tTYPE\\tMOTIF\\tID\\n&/'"

#-----------------------------------------------------------------------------
# call add_dist.py to add 3 columns:
# - distance to next singleton (D1: any individual)
# - distance to next singleton (D2: same individual)
# - number of singletons
#-----------------------------------------------------------------------------
distcmd = config["pythonpath"] + " " + config["sourcedir"] + "scripts/add_dist.py -i - "

#-----------------------------------------------------------------------------
# sort by ID, then by position
#-----------------------------------------------------------------------------
sortcmd = "sort -k7,7 -k2,2n " # + \
    # sourcedir + "topmed_freeze3_singletons/" + \
    # chrom + ".freeze3.singletons.txt"

#-----------------------------------------------------------------------------
# annotate each singleton with 3 columns:
# - cluster ID
# - cluster width
# - number of singletons in cluster
#-----------------------------------------------------------------------------
# clustoutdir = sourcedir + "topmed_freeze3_singletons/" + args.pop + "/sorted/"
# clustoutfile = clustoutdir + chrom + ".freeze3.singletons.sort.txt"
# outfile = sourcedir + "topmed_freeze3_singletons/" + args.pop + "/sorted/" + \
outfile = config["sourcedir"] + "topmed_freeze5_singletons_mask/" + args.pop + "/sorted/" + \
    chrom + ".freeze5.singletons.sort.txt"

# outfile = sourcedir + "topmed_freeze5_dnms/" + args.pop + "/sorted/" + \
#     chrom + ".freeze3.singletons.sort.txt"

clustcmd = config["pythonpath"] + " " + config["sourcedir"] + "process_data/cluster_id.py -i - > " + outfile


#-----------------------------------------------------------------------------
# assemble and run pipe
#-----------------------------------------------------------------------------
cmd = concatcmd + " | " + \
    headercmd + " | " + \
    distcmd + " | " + \
    sortcmd + " | " + \
    clustcmd

# cmd = sortcmd + " | " + clustcmd
eprint(cmd)
os.system(cmd)

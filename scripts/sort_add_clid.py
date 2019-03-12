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

args = parser.parse_args()

###############################################################################
# print to stderr
###############################################################################
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


#-----------------------------------------------------------------------------
# define paths
#-----------------------------------------------------------------------------
chrom = "chr" + str(args.chr)
sourcedir = "/net/snowwhite/home/jedidiah/"
pythonpath = sourcedir + "anaconda3/bin/python "
bcftoolspath = sourcedir + "anaconda3/envs/anno/bin/bcftools "

# vcfdir = "/net/topmed/working/dtaliun/TOPMed_paper_phased/unrelated/"
# vcftail = ".freeze3a.gtonly.Eagle.Phased.unrelated_samples_for_analysis.min_ac_1.vcf.gz"
# vcftest = "~/chr9.freeze3a.test.vcf"

# viewcmd = bcftoolspath + "view -c 1 -C 1 " + vcfdir + chrom + vcftail
# annocmd = pythonpath + sourcedir + "scripts/annotate_motifs.py -i - -f " + sourcedir + "/GRCh37-lite.fa -l 7"
# querycmd = bcftoolspath + "query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/type\\t%INFO/motif\\t%INFO/sample\\n'"
# headercmd = "/bin/sed '1s/.*/CHR\\tPOS\\tREF\\tALT\\tTYPE\\tMOTIF\\tID\\n&/'"

#-----------------------------------------------------------------------------
# concatenate 1Mb chunks from each chromosome
# uses 'ls -v' to list 1Mb chunk files in correct (numeric) order
#-----------------------------------------------------------------------------
# sourcedir + "topmed_freeze3_singletons/" + args.pop + "/" + \
concatcmd = "ls -v " + \
    sourcedir + "topmed_freeze5_singletons_mask/" + args.pop + "/" + \
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
distcmd = pythonpath + sourcedir + "scripts/add_dist.py -i - " 

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
outfile = sourcedir + "topmed_freeze5_singletons_mask/" + args.pop + "/sorted/" + \
    chrom + ".freeze5.singletons.sort.txt"

# outfile = sourcedir + "topmed_freeze5_dnms/" + args.pop + "/sorted/" + \
#     chrom + ".freeze3.singletons.sort.txt"

clustcmd = pythonpath + sourcedir + "scripts/cluster_id.py -i - > " + outfile
    

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

import sys
import os
import argparse
import csv
from collections import defaultdict

###############################################################################
# print to stderr
###############################################################################
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

###############################################################################
# parse args
###############################################################################
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input",
                    help="path to input file (use \"--input -\" \
                        to accept input from stdin)",
                    required=True,
                    nargs='?',
                    type=str,
                    # metavar='',
                    default=sys.stdin)

args = parser.parse_args()

if args.input == "-":
    f = sys.stdin.read().splitlines()
else:
    f = open(args.input, 'r')

###############################################################################
# read input
###############################################################################
# with open(args.input, 'rb') as f:
reader = csv.DictReader(f, delimiter='\t')

prevrow = {'CHR': '', 'ID': '', 'POS': '0'}

site_count = 0
cluster_dict = {}

for row in reader:
    
    # print header
    if site_count == 0:
        print("\t".join(row.keys()) + "\tcl_ID\tcl_LEN\tcl_NUM")
    
    # get site positions for previous and current row
    prevpos = prevrow['POS']
    pos = row['POS']
    
    newrowprev = list(prevrow.values())
    newrow = list(row.values())
    
    # D2 == 0 indicates new ID
    # if site <100bp from previous site in same sample, add to cluster cache
    if int(row['D2']) < 100 and int(row['D2']) > 0:
        
        cluster_dict[prevpos] = newrowprev
        cluster_dict[pos] = newrow
        
        maxpos = max(cluster_dict.keys(), key=(lambda k: cluster_dict[k]))
        minpos = min(cluster_dict.keys(), key=(lambda k: cluster_dict[k]))
        
        cl_ID = str(minpos) + ":" + str(maxpos)
        cl_len = int(maxpos) - int(minpos) + 1
        cl_num = str(len(cluster_dict))

    # if current site >100bp from previous, either print cached cluster
    # or print previous site as unclustered
    else:
        if len(cluster_dict) > 1:
            for key, values in cluster_dict.items():
                newcols = "\t" + cl_ID + "\t" + str(cl_len) + "\t" + cl_num
                print("\t".join(values) + newcols)
        
        # print unclustered with "UC"
        else:
            if int(prevrow['POS']) > 0:
                print("\t".join(newrowprev)+"\tUC\t0\t0")
        
        # reset cluster cache
        cluster_dict = {}
    
    # update prevrow    
    site_count += 1        
    prevrow = row

    # testing
    # if site_count > 200:
    #     break

# ensure last row gets printed    
print("\t".join(newrow)+"\tUC\t0\t0")  

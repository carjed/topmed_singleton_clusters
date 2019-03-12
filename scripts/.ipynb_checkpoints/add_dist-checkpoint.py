import sys
import argparse

import csv
from collections import defaultdict

# outfile = "/mnt/norbert/data/bridges/bridges_ervs_3560_dist.txt"
# outfh = open(outfile, 'w')

# infile = "/mnt/norbert/data/bridges/bridges_ervs_3560.txt"

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input",
                    help="path to input VCF or BCF file (use \"--input -\" \
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

# with open(args.input, 'rb') as f:
reader = csv.DictReader(f, delimiter='\t')

prevrow = {'CHR': '', 'ID': '', 'POS': '0'}
last_site = {}
interval_count = {}

site_count = 0
# D1_count = 0
# int_count = 0
for row in reader:
    ID = row['ID']
    chrom = row['CHR']
    
    # print header and index first item
    if site_count == 0:
        print("\t".join(row.keys()) + "\tD1\tD2\tCOUNT")
        interval_count[ID] = 0
    
    if chrom == prevrow['CHR']:
        D1 = abs(int(row['POS']) - int(prevrow['POS']))
        # D1_count += 1
        if ID in last_site:
            D2 = abs(int(row['POS']) - int(last_site[ID]['POS']))
        else:
            D2 = 0
            interval_count[ID] = 0

        interval_count.update({n: 1 + interval_count[n] for n in interval_count.keys()})
        # if interval_count[ID] 

        if ID != ".":
            int_count = interval_count[ID]-1
            newrow = list(row.values()) + [str(D1), str(D2), str(int_count)]
            print("\t".join(newrow))
            interval_count[ID] = 0

    # ensure first variant per chromosome is printed
    else:
        newrow = list(row.values()) + [str(0), str(0), str(0)]
        print("\t".join(newrow))
        last_site = {}
        interval_count = {}
        interval_count[ID] = 0
        
    site_count += 1        
    prevrow = row
    last_site[ID] = row
    # interval_count[ID] = 
    
        
#         if site_count > 100:
#             break

# outfh.close()

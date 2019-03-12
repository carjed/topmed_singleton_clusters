#!/usr/bin/python

from __future__ import print_function
import os
import sys

sys.path.append(os.getcwd())

import shutil
import textwrap
import argparse
import itertools
import timeit
import time
import multiprocessing
import numpy as np
from joblib import Parallel, delayed
from subprocess import call
from distutils.dir_util import copy_tree

import textwrap
import collections
import csv
import re
from pandas import *
import numpy as np
import cyvcf2 as vcf
from cyvcf2 import VCF, Writer
from scipy.stats import chisquare
from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


###############################################################################
# Parse arguments
###############################################################################
start = timeit.default_timer()

num_cores = multiprocessing.cpu_count()

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input",
                    help="path to input VCF or BCF file (use \"--input -\" \
                        to accept input from stdin)",
                    required=True,
                    nargs='?',
                    type=str,
                    # metavar='',
                    default=sys.stdin)

parser.add_argument("-o", "--output",
                    help="path to output VCF (defaults to stdout)",
                    nargs='?',
                    type=str,
                    # metavar='',
                    default="-")

parser.add_argument("-f", "--fastafile",
                    help="reference fasta file for looking up sequence motifs",
                    nargs='?',
                    type=str,
                    metavar='',
                    default="chr20.fasta.gz")

parser.add_argument("-r", "--ratefile",
                    help="add info field containing relative mutation rate of \
                        the variant",
                    nargs='?',
                    type=str)
                    
parser.add_argument("-s", "--storesample",
                    help="add info field containing ID(s) of sample(s) \
                        carrying the variant",
                    action="store_true")

parser.add_argument("-b", "--bridges",
                    help="fix headers for BRIDGES data",
                    action="store_true")

motif_length_opts = [1,3,5,7]
mlo_str = ",".join(str(x) for x in motif_length_opts)

parser.add_argument("-l", "--length",
                    help="motif length. Allowed values are: " + mlo_str,
                    nargs='?',
                    type=int,
                    choices=motif_length_opts,
                    metavar='',
                    default=7)

parser.add_argument("-v", "--verbose",
                    help="Enable verbose logging",
                    action="store_true")

args = parser.parse_args()

###############################################################################
# print to stderr
###############################################################################
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

###############################################################################
# define k-mer mutation subtypes
###############################################################################
def indexSubtypes(motiflength):
    categories = ["A_C", "A_G", "A_T", "C_G", "C_T", "C_A"]
    bases = ["A", "C", "G", "T"]
    flank = (motiflength-1)//2

    if motiflength > 1:
        kmers = itertools.product(bases, repeat=motiflength-1)

        subtypes_list = []

        for kmer in kmers:
            kmerstr = ''.join(kmer)

            for category in categories:
                ref = category[0]

                subtype = category + "." \
                    + kmerstr[0:flank] + ref + kmerstr[flank:(motiflength-1)]

                subtypes_list.append(subtype)
    else:
        ext = [".A", ".C"]
        extr = list(np.repeat(ext,3))
        subtypes_list = [m+n for m,n in zip(categories,extr)]

    i = 0
    subtypes_dict = {}
    for subtype in sorted(subtypes_list):
        subtypes_dict[subtype] = i
        i += 1

    return subtypes_dict

###############################################################################
# collapse mutation types per strand symmetry
###############################################################################
def getCategory(mu_type):
    if re.match("^[ACGT]*$", mu_type):
        if (mu_type == "AC" or mu_type == "TG"):
            category = "A_C"
        if (mu_type == "AG" or mu_type == "TC"):
            category = "A_G"
        if (mu_type == "AT" or mu_type == "TA"):
            category = "A_T"
        if (mu_type == "CA" or mu_type == "GT"):
            category = "C_A"
        if (mu_type == "CG" or mu_type == "GC"):
            category = "C_G"
        if (mu_type == "CT" or mu_type == "GA"):
            category = "C_T"
    else:
        category = "unknown"
    return category

###############################################################################
# query reference genome for local sequence motif
###############################################################################
def getMotif(pos, sequence):
    motif = Seq(sequence, IUPAC.unambiguous_dna)
    altmotif = motif.reverse_complement()
    central_base = (len(motif)-1)//2

    m1 = motif[central_base]
    m2 = altmotif[central_base]

    if m1 < m2:
        motif_a = motif
    else:
        motif_a = altmotif

    return motif_a

###############################################################################
# build dictionary of rates, indexed by subtype
###############################################################################
def indexRates(ratefile):
    rates_dict =  {}
    with open(ratefile, "r") as f:
        for line in f:
            s = line.strip().split("\t")
            rates_dict[s[0]] = float(s[1])
        return rates_dict

###############################################################################
# Main function for parsing VCF
###############################################################################
def processVCF(args, inputvcf, outputvcf, subtypes_dict):
    if args.verbose:
        eprint("----------------------------------")
        eprint("INITIALIZING REFERENCE GENOME")
        eprint("----------------------------------")
    fasta_reader = Fasta(args.fastafile, read_ahead=1000000)
    eprint("\tDONE") if args.verbose else None
    # record_dict = SeqIO.to_dict(SeqIO.parse(args.fastafile, "fasta"))

    vcf_reader = VCF(inputvcf,
        mode='rb', gts012=True, lazy=True)

    samples = vcf_reader.samples

    vcf_reader.add_info_to_header({'ID': 'type', 
        'Description': 'mutation type',
        'Type':'Character', 
        'Number': '1'})

    vcf_reader.add_info_to_header({'ID': 'motif', 
        'Description': 'k-mer centered at site',
        'Type':'Character', 
        'Number': '1'})

    vcf_reader.add_info_to_header({'ID': 'subtype', 
        'Description': 'k-mer subtype',
        'Type':'Character', 
        'Number': '1'})

    if args.ratefile:
        rates_dict = indexRates(args.ratefile)
        vcf_reader.add_info_to_header({'ID': 'rel_rate', 
            'Description': 'relative mutation rate',
            'Type':'Float', 
            'Number': '1'})

    if args.storesample:
        vcf_reader.add_info_to_header({'ID': 'sample', 
            'Description': 'sample ID',
            'Type':'Character', 
            'Number': '1'})

    # BRIDGES VCFs have several INFO fields that are not defined in the header
    # These lines add the necessary definitions to the output file
    if args.bridges:
        vcf_reader.add_info_to_header({'ID': 'AZ', 'Description': 'AZ',
            'Type':'Character', 'Number': '1'})
        vcf_reader.add_info_to_header({'ID': 'FIC', 'Description': 'FIC',
            'Type':'Character', 'Number': '1'})
        vcf_reader.add_info_to_header({'ID': 'SLRT', 'Description': 'SLRT',
            'Type':'Character', 'Number': '1'})
        vcf_reader.add_info_to_header({'ID': 'LBS', 'Description': 'LBS',
            'Type':'Character', 'Number': '1'})
        vcf_reader.add_info_to_header({'ID': 'OBS', 'Description': 'OBS',
            'Type':'Character', 'Number': '1'})
        vcf_reader.add_info_to_header({'ID': 'LQR', 'Description': 'LQR',
            'Type':'Character', 'Number': '1'})
        vcf_reader.add_info_to_header({'ID': 'SVM', 'Description': 'SVM',
            'Type':'Character', 'Number': '1'})

    # initialize output VCF
    w = Writer(outputvcf, vcf_reader)

    nbp = (args.length-1)//2

    # Query records in VCF and build matrix
    if args.verbose:
        eprint("----------------------------------")
        eprint("PARSING VCF RECORDS")
        eprint("----------------------------------")
    numsites_keep = 0
    numsites_skip = 0
    chrseq = '0'

    for record in vcf_reader:

        # only annotate biallelic SNPs that pass filters
        if record.is_snp and len(record.ALT)==1 and record.FILTER is None:
            # eprint("SNP check: PASS")
            acval = record.INFO['AC']
#             eprint(record.POS, acval)

            # check and update chromosome sequence
            if record.CHROM != chrseq:
                sequence = fasta_reader[record.CHROM]
                chrseq = record.CHROM

            if nbp > 0:
                lseq = sequence[record.POS-(nbp+1):record.POS+nbp].seq
            else:
                lseq = sequence[record.POS-1].seq

            mu_type = record.REF + str(record.ALT[0])
            category = getCategory(mu_type)
            motif_a = getMotif(record.POS, lseq)
            subtype = str(category + "." + motif_a)

            if subtype in subtypes_dict:
                record.INFO["type"] = str(category)

                record.INFO["motif"] = str(motif_a)

                # st = subtypes_dict[subtype]
                record.INFO["subtype"] = subtype
                
                if args.ratefile:
                    record.INFO["rel_rate"] = rates_dict[subtype]

                gt_new = record.gt_types
                if (3 in record.gt_types):
                    gt_new[gt_new == 3] = 0

                # sample = samples[record.gt_types.tolist().index(1)]
                # s1 = ",".join(np.array([samples])(np.nonzero(record.gt_types)[0]).tolist())
                # s1 = record.gt_types.tolist().index(1)
                # eprint(s1)
                # eprint(s1)
                # sample = "test"
                if args.storesample:
                    samples_np = np.array(samples)
                    sample = ",".join(samples_np[np.nonzero(gt_new)[0].tolist()].tolist())
                    record.INFO["sample"] = sample

                # w.write_record(record)
                numsites_keep += 1

            else:
                numsites_skip += 1
                record.INFO["type"] = "."
                record.INFO["motif"] = "."
                record.INFO["subtype"] = "."
                
                if args.storesample:
                    record.INFO["sample"] = "."

            w.write_record(record)

            if args.verbose:
                if (numsites_keep%100000==0):
                    eprint("...", numsites_keep, "sites processed",
                        "(", numsites_skip, "sites skipped)")

    if args.verbose:
        eprint("----------------------------------")
        eprint("VCF PROCESSING COMPLETE")
        eprint("----------------------------------")
        eprint(numsites_keep, "sites kept")
        eprint(numsites_skip, "sites skipped")


    w.close(); vcf_reader.close()

###############################################################################
# index subtypes
###############################################################################
if args.verbose:
    eprint("----------------------------------")
    eprint("INDEXING SUBTYPES")
    eprint("----------------------------------")
subtypes_dict = indexSubtypes(args.length)

if args.verbose:
    eprint("DONE")

###############################################################################
# Annotate VCF
###############################################################################
processVCF(args, args.input, args.output, subtypes_dict)

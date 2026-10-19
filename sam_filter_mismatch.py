#!/usr/bin/env python3

import sys
import argparse
import re

# Argument parsing
parser = argparse.ArgumentParser("Outputs sam lines equal or greater than PID")
parser.add_argument("-a", "--align", dest="align", type=str,  help="Sam file output from bowtie2 of contaminant", nargs="+", required=True)
args = parser.parse_args()

# Generator function for tokenizing file lines
def line_splitter(f, skip=False):
    with open(f) as inputf:
        # If skipping first line
        if skip:
            next(inputf)
        for line in inputf:
            yield line.rstrip().split("\t")

# Calculate percent identity from MD:Z string
def calc_pid(mdz_string, seq_len):
    mat_len = 0

    num_deletes = mdz_string.count("^")
    rem_deletes = mdz_string.replace("^", "")

    list_of_matches = re.split(r"[A,C,T,G,N]", rem_deletes)
    list_of_matches = [x for x in list_of_matches if x]

    for e in list_of_matches:
        mat_len += int(e)
    
    mat_len -= num_deletes

    return mat_len / seq_len

#Parsing file
pid_cutoff = 0.97

for f in args.align:
    for l in line_splitter(f):
        if calc_pid(l[-2].replace("MD:Z:", ""), len(l[9])) >= pid_cutoff:
            print("\t".join(l))

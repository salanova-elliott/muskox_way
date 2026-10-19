#!/usr/bin/env python

import sys
import argparse

# Argument parsing
parser = argparse.ArgumentParser("Computes the proportion of unique k-mers (k=31)")
parser.add_argument("-s", "--sequences", dest="seq", type=str, help="Text file of sequences (column 10 from sam file, unique mappings)", required=True)
args = parser.parse_args()

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

# Generator function for tokenizing file lines
def line_splitter(f):
    with open(f) as inputf:
        for line in inputf:
            yield line.rstrip().split("\t")

def count_kmer(k, seq):
	l = len(seq)
	if l < k: return
	for i in range(l - k + 1):
		kmer_for = seq[i:(i+k)]
		if 'N' in kmer_for: continue
		kmer_rev = kmer_for.translate(comp_tab)[::-1]
		if kmer_for < kmer_rev: kmer = kmer_for
		else: kmer = kmer_rev
		if kmer in counter:
			counter[kmer] += 1
		else: counter[kmer] = 1

def count_file(k, seq_file):
	for line in line_splitter(seq_file):
		count_kmer(k, line[0].upper())
	return

counter = {}
count_file(31, args.seq)

sum_uniq = len(counter)
sum_total = sum(counter.values())

prop_uniq = sum_uniq /sum_total

print(f"{prop_uniq}\t{sum_uniq}\t{sum_total}")
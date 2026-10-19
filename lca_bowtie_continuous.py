#!/usr/bin/env python3

import sys
import argparse
from rocksdict import Rdict, Options, AccessType
import re

# Argument parsing
parser = argparse.ArgumentParser("Outputs LCA for each read input in multimapped SAM format from bowtie2. MUST BE READ SORTED INPUT.")
parser.add_argument("-a", "--align", dest="align", type=str,  help="Sam file output from bowtie2", required=True)
parser.add_argument("-n", "--names", dest="names", type=str, help="NCBI names.dmp file.", required=False, default="/cluster/projects/nn9849k/Reference/taxdump/names.dmp")
parser.add_argument("-s", "--nodes", dest="nodes", type=str, help="NCBI nodes.dmp file.", required=False, default="/cluster/projects/nn9849k/Reference/taxdump/nodes.dmp")
parser.add_argument("-m", "--merged", dest="merged", type=str, help="NCBI merged.dmp file.", required=False, default="/cluster/projects/nn9849k/Reference/taxdump/merged.dmp")
parser.add_argument("-d", "--database", dest="db", type=str, help="Rocksdb of accid to taxid", required=False)
args = parser.parse_args()

# Generator function for tokenizing file lines
def line_splitter(f):
    with open(f) as inputf:
        for line in inputf:
            yield line.rstrip().split("\t")

# Get full taxonomy
def parse_to_root(tid):
    taxonomy = []
    #if tid == "2":
    #    return ["1", "2"]

    while tid != "1":
        taxonomy.append(tid)
        tid = parent_dict[tid][0]
    taxonomy.append("1")
    return taxonomy[::-1]

# Return name at rank
def return_name(tid, taxon_level):
    while tid != "1":
        if parent_dict[tid][1] == taxon_level:
            return name_dict[tid]
        else:
            tid = parent_dict[tid][0]
    return "NA"

# Returns LCA tid from list of lists of taxonomies
def return_lca_tid(t_list):
    min_list_len = len(min(t_list, key=len))
    for i in range(min_list_len):
        if any(e[i] != t_list[0][i] for e in t_list):
            return t_list[0][i-1]
    return t_list[0][min_list_len - 1]

db = Rdict(args.db, access_type=AccessType.read_only())

#Makes a dict of name and taxid {taxid : name}
name_dict = {}
for tokens in line_splitter(args.names):
	#Ignores synonyms, blast names, etc.
	if tokens[6] == "scientific name":
		name_dict[tokens[0]] = tokens[2]
print("Done loading taxid to name data - {}".format(len(name_dict)), file=sys.stderr)

#Makes the parent node dictionary {taxid : (parent_taxid, taxonomic_level)}
parent_dict = {}
for tokens in line_splitter(args.nodes):
	parent_dict[tokens[0]] = (tokens[2], tokens[4])
print("Done loading parent data - {}".format(len(parent_dict)), file=sys.stderr)


#Parsing file
current_read = ""
taxonomy_list = []
tid_set = set()
read_len = 0

#Initialize
for l in line_splitter(args.align):
    current_read = l[0]
    try:
        a_taxid = db.get(bytes(l[2], encoding="utf-8")).decode("utf-8")
    except:
        a_taxid = "-1"
        print(f"NotFound: {l[2]}. Bad initial read", file=sys.stderr)
        break
    tid_set.add(a_taxid)
    taxonomy_list.append(parse_to_root(a_taxid))
    read_len = len(l[9])
    break

# Main loop
r_counter = 0
for l in line_splitter(args.align):

    # If still current read
    if l[0] == current_read:
        try:
            a_taxid = db.get(bytes(l[2], encoding="utf-8")).decode("utf-8")
        except:
            a_taxid = "-1"
            print(f"NotFound: {l[2]}", file=sys.stderr)
            continue
        
        # Catch merged taxids
        try:
            _ = parent_dict[a_taxid]
        except:
			#Looks through merged.dmp for 
            for tokens in line_splitter(args.merged):
                found = False
				#If the taxid is found
                if tokens[0] == a_taxid:
                    a_taxid = tokens[2]
                    found = True
                    break
            if not found:
                print(f"Taxid: {a_taxid} not found in merged.dmp either! Try downloading an updated taxdump", file=sys.stderr)
                break

        if a_taxid not in tid_set:
            taxonomy_list.append(parse_to_root(a_taxid))
            tid_set.add(a_taxid)
        r_counter += 1
    
    # If next read, output current
    else:
        lca_tid = return_lca_tid(taxonomy_list)
            
        print(f"{current_read}\t{lca_tid}\t{read_len}\t{r_counter}", end="")
        for t in parse_to_root(lca_tid)[:0:-1]:
            print(f"\t{parent_dict[t][1]}: {name_dict[t]}", end="")
        print("")

        # Init next read
        try:
            a_taxid = db.get(bytes(l[2], encoding="utf-8")).decode("utf-8")
        except:
            a_taxid = "-1"
            print(f"NotFound: {l[2]}", file=sys.stderr)
            continue

        # Catch merged taxids
        try:
            _ = parent_dict[a_taxid]
        except:
			#Looks through merged.dmp for 
            for tokens in line_splitter(args.merged):
                found = False
				#If the taxid is found
                if tokens[0] == a_taxid:
                    a_taxid = tokens[2]
                    found = True
                    break
            if not found:
                print(f"Taxid: {a_taxid} not found in merged.dmp either! Try downloading an updated taxdump", file=sys.stderr)
                break

        current_read = l[0]
        taxonomy_list = [parse_to_root(a_taxid)]
        r_counter = 1
        tid_set = set(a_taxid)
        read_len = len(l[9])
        lca_tid = 1

#Output final read
lca_tid = return_lca_tid(taxonomy_list)
          
print(f"{current_read}\t{lca_tid}\t{read_len}\t{r_counter}", end="")
for t in parse_to_root(lca_tid)[:0:-1]:
    print(f"\t{parent_dict[t][1]}: {name_dict[t]}", end="")

#!/usr/bin/env python3
import re
import sys
import random
import os
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import yaml
from pprint import pprint

# Ensure logs directory exists
os.makedirs("logs", exist_ok=True)
log_file = "logs/dictgen.log"

# Read configuration
try:
    with open("config.yaml", "r") as f:
        config = yaml.safe_load(f)
except Exception as e:
    with open(log_file, "w") as log:
        log.write(f"Error reading config.yaml: {str(e)}\n")
    sys.exit(f"Error reading config.yaml: {str(e)}")

tree_topology = config.get("tree_topology")
if not tree_topology:
    with open(log_file, "w") as log:
        log.write("Error: 'tree_topology' must be specified in config.yaml\n")
    sys.exit("Error: 'tree_topology' must be specified in config.yaml")

# Function to extract species from Newick tree
def extract_species(newick_tree):
    newick_tree = re.sub(r':[\d.]+', '', newick_tree)
    newick_tree = newick_tree.replace(' ', '')
    species = re.findall(r'[\w]+', newick_tree)
    return sorted(list(set(species)))

# Get all species
all_species = extract_species(tree_topology)

# Log configuration
with open(log_file, "w") as log:
    log.write(f"Tree topology: {tree_topology}\n")
    log.write(f"All species: {all_species}\n")

# Ensure output directory
os.makedirs("ref", exist_ok=True)

# Create empty reference.list
with open("reference.list", "w") as f:
    f.write("")

# Read input files
keys = "goodalignments.txt"
values = sys.argv[1]

try:
    with open(keys) as f:
        keylist = [line.strip() for line in f if line.strip()]
    with open(values) as f:
        valuelist = [line.strip() for line in f if line.strip()]
except Exception as e:
    with open(log_file, "a") as log:
        log.write(f"Error reading input files ({keys}, {values}): {str(e)}\n")
    sys.exit(f"Error reading input files: {str(e)}")

# Validate and extract basenames from keylist
basenames = []
for key in keylist:
    filename = os.path.basename(key)
    if not filename.endswith(".fa.prunned"):
        with open(log_file, "a") as log:
            log.write(f"Warning: Invalid filename format in {key}, expected .fa.prunned\n")
        continue
    base = filename.replace(".fa.prunned", "")
    basenames.append(base)

# Log inputs
with open(log_file, "a") as log:
    log.write(f"Keylist: {keylist}\n")
    log.write(f"Basenames: {basenames}\n")
    log.write(f"Valuelist: {valuelist}\n")

# Helper functions
def resolve_key_path(filename):
    if os.path.exists(filename):
        return filename
    chrom = filename.split('.')[0]
    possible = os.path.join("query", chrom, filename)
    if os.path.exists(possible):
        return possible
    with open(log_file, "a") as log:
        log.write(f"Warning: Key file {filename} not found\n")
    return filename

def resolve_neutral_path(filename):
    if os.path.exists(filename):
        return filename
    possible = os.path.join("neutral_set", filename)
    if os.path.exists(possible):
        return possible
    with open(log_file, "a") as log:
        log.write(f"Warning: Neutral file {filename} not found\n")
    return filename

# Process replicates
for replicate in range(10):
    dictionary = {}
    for key in keylist:
        random.shuffle(valuelist)
        dictionary[key] = valuelist[:20]
    
    with open("global.dict", "w") as fielddict_file:
        pprint(dictionary, fielddict_file)
    with open(log_file, "a") as log:
        log.write(f"Replicate {replicate}: Wrote global.dict\n")
    
    for key in keylist:
        filename = os.path.basename(key)
        if not filename.endswith(".fa.prunned"):
            continue
        base = filename.replace(".fa.prunned", "")
        neutral_list = dictionary[key]
        n = 0
        combined_seq = MultipleSeqAlignment([
            SeqRecord(Seq(''), id=species) for species in all_species
        ])
        combined_seq.sort()
        
        for ref in neutral_list:
            n += 1
            aligned_file = resolve_neutral_path(ref)
            try:
                seq_records = AlignIO.read(aligned_file, 'fasta')
                seq_species = sorted([rec.id for rec in seq_records])
                if seq_species != all_species:
                    with open(log_file, "a") as log:
                        log.write(f"Warning: Species mismatch in {aligned_file}: expected {all_species}, got {seq_species}\n")
                    continue
                for rec in seq_records:
                    rec.description = rec.id
                seq_records.sort()
                combined_seq = combined_seq + seq_records
                combined_seq.description = ""
            except Exception as e:
                with open(log_file, "a") as log:
                    log.write(f"Error processing {aligned_file}: {str(e)}\n")
                continue
        
        for rec in combined_seq:
            rec.description = rec.id
        
        out_filename = os.path.join("ref", f"{base}.{replicate}.ref")
        try:
            with open(out_filename, 'w') as write_file:
                AlignIO.write(combined_seq, write_file, 'fasta')
            with open(log_file, "a") as log:
                log.write(f"Replicate {replicate}: Wrote {out_filename}\n")
        except Exception as e:
            with open(log_file, "a") as log:
                log.write(f"Error writing {out_filename}: {str(e)}\n")
            continue
        
        with open('reference.list', 'a') as referencelist:
            referencelist.write(f"{key}\t{n}\t{replicate}\n")

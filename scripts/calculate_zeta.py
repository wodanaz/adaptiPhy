#!/usr/bin/env python3
import re
import csv
import os
import fnmatch
import argparse
import yaml
import sys
from dendropy import Tree
# ------------------------------------------------------------------
# Force correct working directory (project root)
# ------------------------------------------------------------------
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(script_dir, ".."))
os.chdir(project_root)
# ------------------------------------------------------------------
# Load config
# ------------------------------------------------------------------
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)
foreground_branches = config.get("foreground_branches", [])
internal_branches = config.get("internal_branches", [])  # List of dicts: {'name': str, 'tips': list[str]}
if not foreground_branches:
    sys.exit("Error: 'foreground_branches' must be specified and non-empty in config.yaml")
# ------------------------------------------------------------------
# Function to parse TREE line and extract branch lengths
# ------------------------------------------------------------------
def parse_tree(tree_line):
    branch_lengths = {}
    # Match species:branch_length (e.g., hg19:0.00806296)
    matches = re.findall(r'([A-Za-z0-9_]+):([0-9.eE+-]+)', tree_line)
    for branch, length in matches:
        try:
            branch_lengths[branch] = float(length)
        except ValueError:
            continue
    return branch_lengths
# ------------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Calculate zeta values from phyloFit .mod files")
parser.add_argument("ref_list", help="Path to reference.list")
parser.add_argument("query_dir", help="Directory with query .mod files")
parser.add_argument("ref_dir", help="Directory with reference .mod files")
parser.add_argument("summary_table", help="Output summary table path")
parser.add_argument("--log", default="logs/calculate_zeta.log", help="Log file")
parser.add_argument("--chrom", default=None, help="Optional chromosome to filter basenames (e.g., 'chr1')")
args = parser.parse_args()
ref_list = args.ref_list
query_dir = args.query_dir
ref_dir = args.ref_dir
summary_table = args.summary_table
log_file = args.log
chrom_filter = args.chrom
os.makedirs(os.path.dirname(log_file), exist_ok=True)
with open(log_file, "w") as log:
    log.write("=== calculate_zeta.py started ===\n")
    log.write(f"Query dir: {query_dir}\n")
    log.write(f"Ref dir: {ref_dir}\n")
    log.write(f"Foreground branches: {foreground_branches}\n")
    if chrom_filter:
        log.write(f"Filtering for chromosome: {chrom_filter}\n")
# ------------------------------------------------------------------
# Extract unique basenames from reference.list (column 1 = full path or basename)
# ------------------------------------------------------------------
bases = []
seen = set()
with open(ref_list, "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 1:
            continue
        full_path = parts[0]
        base = os.path.basename(full_path).replace(".fa", "")
        if base not in seen:
            seen.add(base)
            bases.append(base)
# Filter by chrom if specified
if chrom_filter:
    bases = [b for b in bases if b.startswith(f"{chrom_filter}.")]
with open(log_file, "a") as log:
    log.write(f"Found {len(bases)} unique regions\n")
if not bases:
    with open(log_file, "a") as log:
        log.write("No bases found (or after filtering) — created empty output with header\n")
    with open(summary_table, "w", newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["NAME", "BRANCH", "REPL", "zeta"], delimiter="\t")
        writer.writeheader()
else:
    # ------------------------------------------------------------------
    # Pre-load all ref .mod paths into a dict {base: {replicate: path}}
    # ------------------------------------------------------------------
    ref_mod_dict = {}
    for filename in os.listdir(ref_dir):
        if filename.endswith(".mod"):
            parts = filename.split(".")
            if len(parts) >= 3:
                base = ".".join(parts[:-2])  # Everything before .replicate.mod
                replicate = parts[-2]
                path = os.path.join(ref_dir, filename)
                if base not in ref_mod_dict:
                    ref_mod_dict[base] = {}
                ref_mod_dict[base][replicate] = path
    with open(log_file, "a") as log:
        log.write(f"Pre-loaded {len(ref_mod_dict)} unique bases from ref_dir\n")
    # ------------------------------------------------------------------
    # Prepare CSV writer for incremental writing
    # ------------------------------------------------------------------
    headers = ["NAME", "BRANCH", "REPL", "zeta"]
    with open(summary_table, "w", newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers, delimiter="\t")
        writer.writeheader()
    row_count = 0
    # ------------------------------------------------------------------
    # Process each base incrementally
    # ------------------------------------------------------------------
    for base in bases:
        query_mod = os.path.join(query_dir, f"{base}.mod")
        ref_mods = ref_mod_dict.get(base, {})
        if not os.path.exists(query_mod):
            with open(log_file, "a") as log:
                log.write(f"Missing query model: {query_mod}\n")
            continue
        if not ref_mods:
            with open(log_file, "a") as log:
                log.write(f"No reference models found for {base}\n")
            continue
        # Parse query TREE
        query_bl = {}
        query_tree = None
        try:
            with open(query_mod) as f:
                for line in f:
                    if line.startswith("TREE:"):
                        newick = line.split("TREE: ", 1)[1].strip()
                        if not newick.endswith(";"):
                            newick += ";"
                        query_tree = Tree.get_from_string(newick, schema="newick")
                        query_bl = parse_tree(line.strip())
                        break
        except Exception as e:
            with open(log_file, "a") as log:
                log.write(f"Failed to read query {query_mod}: {e}\n")
            continue
        if not query_bl:
            with open(log_file, "a") as log:
                log.write(f"No TREE line found in query {query_mod}\n")
            continue
        # Parse each reference replicate
        for replicate, ref_mod in ref_mods.items():
            try:
                ref_bl = {}
                ref_tree = None
                with open(ref_mod) as f:
                    for line in f:
                        if line.startswith("TREE:"):
                            newick = line.split("TREE: ", 1)[1].strip()
                            if not newick.endswith(";"):
                                newick += ";"
                            ref_tree = Tree.get_from_string(newick, schema="newick")
                            ref_bl = parse_tree(line.strip())
                            break
                    else:
                        continue
            except Exception as e:
                with open(log_file, "a") as log:
                    log.write(f"Failed to read ref {ref_mod}: {e}\n")
                continue
            if not ref_bl:
                with open(log_file, "a") as log:
                    log.write(f"No TREE line found in ref {ref_mod}\n")
                continue
            # Calculate zeta for each foreground branch and write row immediately
            with open(summary_table, "a", newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=headers, delimiter="\t")
                for branch in foreground_branches:
                    q_len = query_bl.get(branch)
                    r_len = ref_bl.get(branch)
                    if q_len is not None and r_len is not None and r_len > 0:
                        zeta = q_len / r_len
                        row = {
                            "NAME": base,
                            "BRANCH": branch,
                            "REPL": replicate,
                            "zeta": zeta
                        }
                        writer.writerow(row)
                        row_count += 1
                    else:
                        with open(log_file, "a") as log:
                            log.write(f"Missing branch length for {branch} in {base} (query: {q_len}, ref: {r_len})\n")
            # Calculate zeta for internal branches (e.g., LCA of hg19 and panTro4)
            for int_branch in internal_branches:
                int_name = int_branch['name']
                tips = int_branch['tips']
                try:
                    # Query tree LCA branch length
                    q_lca = query_tree.mrca(taxon_labels=tips)
                    q_len = q_lca.edge_length if q_lca.edge_length is not None else None
                    # Ref tree LCA branch length
                    r_lca = ref_tree.mrca(taxon_labels=tips)
                    r_len = r_lca.edge_length if r_lca.edge_length is not None else None
                    if q_len is not None and r_len is not None and r_len > 0:
                        zeta = q_len / r_len
                        row = {
                            "NAME": base,
                            "BRANCH": int_name,
                            "REPL": replicate,
                            "zeta": zeta
                        }
                        writer.writerow(row)
                        row_count += 1
                    else:
                        with open(log_file, "a") as log:
                            log.write(f"Missing internal branch length for {int_name} in {base} (query: {q_len}, ref: {r_len})\n")
                except Exception as e:
                    with open(log_file, "a") as log:
                        log.write(f"Failed to compute LCA for {int_name} in {base}: {e}\n")
# ------------------------------------------------------------------
# Final log
# ------------------------------------------------------------------
with open(log_file, "a") as log:
    log.write(f"Wrote {row_count} zeta values to {summary_table}\n")
    log.write("=== calculate_zeta.py finished successfully ===\n")
print(f"calculate_zeta complete — {row_count} rows written to {summary_table}")

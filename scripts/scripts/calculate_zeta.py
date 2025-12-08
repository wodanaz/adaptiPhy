#!/usr/bin/env python3
import re
import pandas as pd
import os
import glob
import argparse
import yaml
import sys

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
args = parser.parse_args()

ref_list = args.ref_list
query_dir = args.query_dir
ref_dir = args.ref_dir
summary_table = args.summary_table
log_file = args.log

os.makedirs(os.path.dirname(log_file), exist_ok=True)

with open(log_file, "w") as log:
    log.write("=== calculate_zeta.py started ===\n")
    log.write(f"Query dir: {query_dir}\n")
    log.write(f"Ref dir: {ref_dir}\n")
    log.write(f"Foreground branches: {foreground_branches}\n")

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

with open(log_file, "a") as log:
    log.write(f"Found {len(bases)} unique regions: {', '.join(bases)}\n")

if not bases:
    sys.exit("No bases found in reference.list")

# ------------------------------------------------------------------
# Process each base
# ------------------------------------------------------------------
zeta_data = []

for base in bases:
    query_mod = os.path.join(query_dir, f"{base}.mod")
    ref_mod_pattern = os.path.join(ref_dir, f"{base}.*.mod")
    ref_mods = sorted(glob.glob(ref_mod_pattern))

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
    try:
        with open(query_mod) as f:
            for line in f:
                if line.startswith("TREE:"):
                    query_bl = parse_tree(line.strip())
                    break
    except Exception as e:
        with open(log_file, "a") as log:
            log.write(f"Failed to read query {query_mod}: {e}\n")
        continue

    # Parse each reference replicate
    for ref_mod in ref_mods:
        try:
            replicate = os.path.basename(ref_mod).split(".")[-2]  # e.g., "8" from ...8.mod
            with open(ref_mod) as f:
                for line in f:
                    if line.startswith("TREE:"):
                        ref_bl = parse_tree(line.strip())
                        break
                else:
                    continue
        except Exception as e:
            with open(log_file, "a") as log:
                log.write(f"Failed to read ref {ref_mod}: {e}\n")
            continue

        # Calculate zeta for each foreground branch
        for branch in foreground_branches:
            q_len = query_bl.get(branch)
            r_len = ref_bl.get(branch)
            if q_len is not None and r_len is not None and r_len > 0:
                zeta = q_len / r_len
            else:
                zeta = "NA"

            zeta_data.append({
                "NAME": base,
                "BRANCH": branch,
                "REPL": replicate,
                "zeta": zeta
            })

# ------------------------------------------------------------------
# Write output
# ------------------------------------------------------------------
df = pd.DataFrame(zeta_data)
df.to_csv(summary_table, sep="\t", index=False)

with open(log_file, "a") as log:
    log.write(f"Wrote {len(df)} zeta values to {summary_table}\n")
    log.write("=== calculate_zeta.py finished successfully ===\n")

print(f"calculate_zeta complete — {len(df)} rows written to {summary_table}")

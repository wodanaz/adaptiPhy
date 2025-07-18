#!/usr/bin/env python3
import re
import pandas as pd
import os
import glob
import argparse
import yaml
import sys

# Read configuration
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)
foreground_branches = config["foreground_branches"]

# Validate that foreground_branches is not empty
if not foreground_branches:
    sys.exit("Error: 'foreground_branches' must be specified in config.yaml")

# Function to parse TREE line and extract branch lengths
def parse_tree(tree_line):
    branch_lengths = {}
    matches = re.findall(r'(\w+):([\d.]+)', tree_line)
    for branch, length in matches:
        branch_lengths[branch] = float(length)
    return branch_lengths

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Calculate zeta values for all bases from query and reference .mod files.")
parser.add_argument("ref_list", help="Path to ref.list file containing base names")
parser.add_argument("query_dir", help="Directory containing query .mod files (e.g., MODELS_HKY85/query)")
parser.add_argument("ref_dir", help="Directory containing reference .mod files (e.g., MODELS_HKY85/ref)")
parser.add_argument("summary_table", help="Path to the input/output summary table")
parser.add_argument("--log", default="calculate_zeta.log", help="Path to the log file")
args = parser.parse_args()

# Assign inputs
ref_list = args.ref_list
query_dir = args.query_dir
ref_dir = args.ref_dir
summary_table = args.summary_table
log_file = args.log

# Read base names from ref.list
with open(ref_list, 'r') as f:
    bases = [line.strip() for line in f if line.strip()]

# Log input details
with open(log_file, 'w') as log:
    log.write(f"Reference list: {ref_list}\n")
    log.write(f"Query directory: {query_dir}\n")
    log.write(f"Reference directory: {ref_dir}\n")
    log.write(f"Found {len(bases)} bases: {', '.join(bases)}\n")
    log.write(f"Foreground branches: {', '.join(foreground_branches)}\n")

# Initialize zeta data
zeta_data = []

# Process each base
for base in bases:
    query_mod = os.path.join(query_dir, f"{base}.mod")
    ref_mods = sorted(glob.glob(os.path.join(ref_dir, f"{base}.*.mod")))  # Match all replicates

    # Validate inputs
    if not os.path.exists(query_mod):
        with open(log_file, 'a') as log:
            log.write(f"Error: Query file {query_mod} does not exist for base {base}\n")
        continue

    if not ref_mods:
        with open(log_file, 'a') as log:
            log.write(f"Error: No reference files found for base {base} in {ref_dir}\n")
        continue

    with open(log_file, 'a') as log:
        log.write(f"Processing base {base} with {len(ref_mods)} reference files:\n")
        for f in ref_mods:
            log.write(f"  {f}\n")

    # Parse query .mod file
    try:
        with open(query_mod, 'r') as f:
            for line in f:
                if line.startswith('TREE:'):
                    query_branch_lengths = parse_tree(line.strip())
                    with open(log_file, 'a') as log:
                        log.write(f"Query TREE for {base}: {line.strip()}\n")
                    break
            else:
                with open(log_file, 'a') as log:
                    log.write(f"Error: No TREE line found in {query_mod}\n")
                continue
    except Exception as e:
        with open(log_file, 'a') as log:
            log.write(f"Error: Failed to read {query_mod}: {str(e)}\n")
        continue

    # Parse reference .mod files and compute zeta per replicate
    valid_ref_files = 0
    for ref_mod in ref_mods:
        try:
            with open(ref_mod, 'r') as f:
                for line in f:
                    if line.startswith('TREE:'):
                        ref_branch_lengths = parse_tree(line.strip())
                        # Extract replicate number from filename (e.g., '0' from 'chr9.53180570-53181126.0.mod')
                        replicate = os.path.basename(ref_mod).split('.')[-2]
                        with open(log_file, 'a') as log:
                            log.write(f"Processed {ref_mod} TREE: {line.strip()}\n")
                        valid_ref_files += 1
                        # Calculate zeta for each branch
                        for branch in foreground_branches:
                            if branch in query_branch_lengths and branch in ref_branch_lengths:
                                zeta = query_branch_lengths[branch] / ref_branch_lengths[branch]
                                zeta_data.append({
                                    'NAME': base,
                                    'REPL': replicate,
                                    'BRANCH': branch,
                                    'zeta': zeta
                                })
                                with open(log_file, 'a') as log:
                                    log.write(f"Zeta {branch} for {base} replicate {replicate}: {query_branch_lengths[branch]} / {ref_branch_lengths[branch]} = {zeta}\n")
                            else:
                                zeta_data.append({
                                    'NAME': base,
                                    'REPL': replicate,
                                    'BRANCH': branch,
                                    'zeta': 'NA'
                                })
                                with open(log_file, 'a') as log:
                                    log.write(f"Warning: Missing branch {branch} in {ref_mod} or query for {base}\n")
                        break
                else:
                    with open(log_file, 'a') as log:
                        log.write(f"Warning: No TREE line found in {ref_mod}\n")
        except Exception as e:
            with open(log_file, 'a') as log:
                log.write(f"Warning: Failed to process {ref_mod}: {str(e)}\n")

    if valid_ref_files == 0:
        with open(log_file, 'a') as log:
            log.write(f"Error: No valid reference files with TREE line found for {base}\n")

# Read or create summary table
if os.path.exists(summary_table) and os.path.getsize(summary_table) > 0:
    try:
        summary_df = pd.read_csv(summary_table, sep='\t')
        with open(log_file, 'a') as log:
            log.write(f"Summary table columns: {summary_df.columns.tolist()}\n")
        if not {'NAME', 'BRANCH', 'REPL', 'zeta'}.issubset(summary_df.columns):
            with open(log_file, 'a') as log:
                log.write(f"Warning: Required columns missing in {summary_table}. Creating new DataFrame.\n")
            summary_df = pd.DataFrame(columns=['NAME', 'BRANCH', 'REPL', 'zeta'])
    except (pd.errors.EmptyDataError, pd.errors.ParserError) as e:
        with open(log_file, 'a') as log:
            log.write(f"Warning: Failed to parse {summary_table}: {str(e)}. Creating new DataFrame.\n")
        summary_df = pd.DataFrame(columns=['NAME', 'BRANCH',  'REPL', 'zeta'])
else:
    with open(log_file, 'a') as log:
        log.write(f"Info: {summary_table} does not exist or is empty. Creating new DataFrame.\n")
    summary_df = pd.DataFrame(columns=['NAME', 'BRANCH',  'REPL',  'zeta'])

# Create new DataFrame from zeta data
new_df = pd.DataFrame(zeta_data)

# Merge with existing summary table
for _, row in new_df.iterrows():
    mask = (summary_df['NAME'] == row['NAME']) & (summary_df['BRANCH'] == row['BRANCH']) & (summary_df['REPL'] == row['REPL'])
    if mask.any():
        summary_df.loc[mask, 'zeta'] = row['zeta']
    else:
        summary_df = pd.concat([summary_df, pd.DataFrame([row])], ignore_index=True)

# Save updated summary table
summary_df.to_csv(summary_table, sep='\t', index=False)

# Log completion
with open(log_file, 'a') as log:
    log.write(f"Updated summary table with zeta values for {len(bases)} bases, total {len(zeta_data)} rows\n")

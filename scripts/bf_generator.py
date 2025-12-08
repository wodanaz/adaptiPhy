#!/usr/bin/env python3
import sys
import random
import os
import yaml

# 1. ALWAYS run from the project root (where config.yaml and Snakefile live)
script_dir = os.path.dirname(os.path.abspath(__file__))      # .../scripts
project_root = os.path.abspath(os.path.join(script_dir, ".."))  # one level up
os.chdir(project_root)

if len(sys.argv) < 2:
    print("Usage: python scripts/bf_generator.py <reference.list>")
    sys.exit(1)

ref_list_file = sys.argv[1]

# 2. Load config safely
config_path = os.path.join(project_root, "config.yaml")
if not os.path.exists(config_path):
    sys.exit(f"ERROR: config.yaml not found at {config_path}")

with open(config_path) as f:
    config = yaml.safe_load(f)

tree_topology = config.get("tree_topology")
foreground_branches = config.get("foreground_branches", [])

if not tree_topology:
    sys.exit("ERROR: 'tree_topology' missing in config.yaml")
if not foreground_branches:
    sys.exit("ERROR: 'foreground_branches' empty in config.yaml")

# 3. Parse reference.list — column 1 can be full path or clean name
entries = []
with open(ref_list_file) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        col1 = parts[0]
        base = os.path.basename(col1).replace(".fa", "")   # works for both full path and clean name
        replicate = parts[2]
        entries.append((base, replicate))

if not entries:
    sys.exit(f"ERROR: No valid entries found in {ref_list_file}")

# 4. Generate .bf files
output_dir = "bf_dir"
os.makedirs(output_dir, exist_ok=True)

models = ["null", "alt"]
total = 0

for base, rep in entries:
    for model in models:
        for branch in foreground_branches:
            seed = random.randint(1, 999999)
            bf_name = f"{base}.{branch}.{model}.{rep}.bf"
            bf_path = os.path.join(output_dir, bf_name)

            with open(bf_path, "w") as bf:
                bf.write(f"random_seed = {seed};\n")
                bf.write(f'quer_seq_file = "good_alignments/{base}.fa";\n')
                bf.write(f'ref_seq_file = "ref/{base}.{rep}.ref";\n')
                bf.write("fit_repl_count = 20;\n")
                bf.write(f'tree = "{tree_topology}";\n')
                bf.write(f'fgrnd_branch_name = "{branch}";\n')
                bf.write(f'res_file = "../res/{base}.{branch}.{model}.{rep}.res";\n')
                bf.write(f'#include "{model}4-fgrnd_spec.model";\n')
            total += 1

print(f"bf_generator.py: Successfully created {total} .bf files in {os.path.abspath(output_dir)}")

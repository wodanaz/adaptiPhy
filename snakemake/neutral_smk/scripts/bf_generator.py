# /cwork/ab620/adaptify_snake/scripts/bf_generator.py
#!/usr/bin/env python3
import sys
import random
import os
import yaml

if len(sys.argv) < 2:
    sys.exit("Usage: python bf_generator.py <ref_list_file>")

# Read configuration
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)
tree_topology = config["tree_topology"]
foreground_branches = config["foreground_branches"]

# Validate that foreground_branches is not empty
if not foreground_branches:
    sys.exit("Error: 'foreground_branches' must be specified in config.yaml")

# Read base names
with open(sys.argv[1]) as f:
    base_names = [line.strip() for line in f if line.strip()]

models = ["null", "alt"]
num_replicates = 10
output_dir = "bf_dir"

os.makedirs(output_dir, exist_ok=True)

for replicate in range(num_replicates):
    for model in models:
        for branch in foreground_branches:
            for base in base_names:
                chrom = base.split('.')[0]
                seed = random.randint(1, 1000)
                out_filename = f"{base}.{branch}.{model}.{replicate}.bf"
                out_path = os.path.join(output_dir, out_filename)
                with open(out_path, "w") as out_file:
                    out_file.write(f"random_seed={seed};\n")
                    out_file.write(f'quer_seq_file = "query/{base}.fa.prunned";\n')
                    out_file.write(f'ref_seq_file = "ref/{base}.{replicate}.ref";\n')
                    out_file.write("fit_repl_count = 20;\n")
                    out_file.write(f'tree = "{tree_topology}";\n')
                    out_file.write(f'fgrnd_branch_name = "{branch}";\n')
                    out_file.write(f'res_file = "../res/{base}.{branch}.{model}.{replicate}.res";\n')
                    out_file.write(f'#include "{model}4-fgrnd_spec.model";\n')

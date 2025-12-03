#!/usr/bin/env python3
"""
scripts/parse_neutral.py

Parse PhyloFit output into filtered tables.

Features:
- Reads multiple tree topologies and foreground branches from config.yaml (accepts Newick with or without end semicolon).
- Parses lines like: MODELS_HKY85/batch*/<region>.mod:TREE: <newick>
- Extracts tip branch lengths, internal edge lengths, total tree length, and computes tip relative branch lengths.
- Writes a permanent full (unfiltered) table: neutral_table_full.tsv.
- Applies tip filter (all tip lengths > 0.001) and IQR filter on foreground relative branch length to produce the final neutral_table and neutralset list.
- Logs progress and sample parse failures to a log file.
"""
import os
import re
import sys
import time
import argparse
import traceback
import yaml
import pandas as pd
from dendropy import Tree

##########################################################################

def load_config():
    """Load config.yaml (must be in cwd) and validate tree_topology/foreground_branch.
    Accepts tree_topology with or without a trailing semicolon and tries again
    if initial parsing fails.
    """
    cfg_path = "config.yaml"
    if not os.path.exists(cfg_path):
        raise FileNotFoundError(f"{cfg_path} not found in the current working directory.")
    with open(cfg_path, "r") as f:
        config = yaml.safe_load(f) or {}

    tree_topology = config.get(
        "tree_topology",
        "(E,(D,(C,(A,B))));"
    )
    foreground_branch = config.get("foreground_branch", "A")

    if not isinstance(tree_topology, str):
        raise ValueError("tree_topology in config.yaml must be a Newick string.")

    try:
        topo_tree = Tree.get_from_string(tree_topology, schema="newick")
    except Exception as first_err:
        if not tree_topology.strip().endswith(";"):
            try:
                tree_topology = tree_topology.strip() + ";"
                topo_tree = Tree.get_from_string(tree_topology, schema="newick")
            except Exception as e2:
                raise ValueError(f"Invalid tree_topology in config.yaml: {e2}")
        else:
            raise ValueError(f"Invalid tree_topology in config.yaml: {first_err}")

    branches = []
    for node in topo_tree.leaf_nodes():
        if node.taxon and node.taxon.label:
            branches.append(node.taxon.label)
    if not branches:
        raise ValueError("No tip labels found in tree_topology.")
    if foreground_branch not in branches:
        print(f"Warning: foreground_branch '{foreground_branch}' not in tree. Using '{branches[0]}'")
        foreground_branch = branches[0]

    return branches, foreground_branch, tree_topology


def log_message(log_file, msg, mode="a"):
    """Append timestamped message to log file (creates directory if needed)."""
    d = os.path.dirname(log_file) or "."
    os.makedirs(d, exist_ok=True)
    with open(log_file, mode) as f:
        f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n")
        f.flush()


def get_internal_edges(tree):
    """Collect internal edge lengths with descriptive labels."""
    internals = {}
    i = 1
    for node in tree.internal_nodes():
        if node.parent_node and node.edge_length is not None:
            children = [n.taxon.label for n in node.child_nodes()
                        if n.is_leaf() and n.taxon and n.taxon.label]
            if children:
                label = f"node_{i}_clade_{'-'.join(sorted(children))}"
            else:
                label = f"node_{i}"
            internals[label] = node.edge_length
            i += 1
    return internals


def parse_line(line, branches):
    """Parse one line from output.hky85.neutral.txt → dict of values.
    Returns None on parse failure or if values are unreliable.
    """
    line = line.strip()
    if not line or not line.startswith("MODELS_HKY85/") or ":TREE:" not in line:
        return None

    parts = line.split(":TREE:", 1)
    if len(parts) != 2:
        return None

    region = parts[0].split("/")[-1].replace(".mod", "")
    newick = parts[1].strip()

    if not newick.endswith(";"):
        newick = newick + ";"

    try:
        t = Tree.get_from_string(newick, schema="newick")
    except Exception:
        return None

    branch_lengths = {}
    try:
        for br in branches:
            node = t.find_node_with_taxon_label(br)
            if node is None:
                return None
            length = node.edge_length
            if length is None or (isinstance(length, float) and length < 1e-10):
                return None
            branch_lengths[br] = length
    except Exception:
        return None

    internal_lengths = get_internal_edges(t)

    total = 0.0
    try:
        for e in t.edges():
            le = getattr(e, "length", None)
            if le is not None and isinstance(le, (int, float)) and le > 1e-10:
                total += le
    except Exception:
        return None

    if total <= 0.0:
        return None

    rel_branches = {f"{b}_rb": (v / total) for b, v in branch_lengths.items()}

    rec = {
        "region": region,
        "total": total,
        **branch_lengths,
        **internal_lengths,
        **rel_branches,
    }
    return rec


##################################################################

def main():
    parser = argparse.ArgumentParser(
        description="Parse neutral tree data into filtered tables."
    )
    parser.add_argument("input_txt", help="e.g. output.hky85.neutral.txt")
    parser.add_argument("output_table", help="neutral_table.txt")
    parser.add_argument("neutral_list", help="neutralset.txt")
    parser.add_argument("--log", default="logs/parse_neutral.log", help="Log file")
    parser.add_argument("--full-table", default=None, help="Optional full (unfiltered) table path; if omitted, defaults to <output_table>_full.tsv")
    args = parser.parse_args()

    log_file     = os.path.abspath(args.log)
    input_txt    = os.path.abspath(args.input_txt)
    output_table = os.path.abspath(args.output_table)
    neutral_list = os.path.abspath(args.neutral_list)
    full_table   = os.path.abspath(args.full_table) if args.full_table else None

    try:
        branches, fg_branch, topo = load_config()
        log_message(log_file, "Starting parse_neutral.py", "w")
        log_message(log_file, f"Input file: {input_txt}")
        log_message(log_file, f"Branches: {', '.join(branches)} | Foreground: {fg_branch}")

        data = []
        skipped = 0
        max_show_skips = 5
        with open(input_txt, "r") as f:
            for line_num, line in enumerate(f, 1):
                parsed = parse_line(line, branches)
                if parsed:
                    data.append(parsed)
                else:
                    skipped += 1
                    if skipped <= max_show_skips:
                        log_message(log_file, f"Skipped line {line_num}: {line.strip()[:200]}")

        log_message(log_file, f"Parsed {len(data)} regions, skipped {skipped}")

        if not data:
            log_message(log_file, "No valid trees – creating empty outputs.")
            open(output_table, "w").close()
            open(neutral_list, "w").close()
            # write empty full (unfiltered) table as permanent file
            default_full = os.path.splitext(output_table)[0] + "_full.tsv"
            full_path = full_table if full_table else default_full
            open(full_path, "w").close()
            log_message(log_file, f"Wrote empty full table → {full_path}")
            return

        df = pd.DataFrame(data).round(8)

        # ---- Save full (unfiltered) table ----
        default_full = os.path.splitext(output_table)[0] + "_full.tsv"
        full_path = full_table if full_table else default_full
        try:
            df.to_csv(full_path, sep="\t", index=False)
            log_message(log_file, f"Full (unfiltered) table saved ({len(df)} rows) → {full_path}")
        except Exception as e:
            log_message(log_file, f"WARNING: could not write full table {full_path}: {e}")

        # ---- Tip filter (>0.001) ----
        tip_cols = branches
        missing = [c for c in tip_cols if c not in df.columns]
        if missing:
            raise ValueError(f"Missing expected tip columns in parsed table: {missing}")

        filtered = df[(df[tip_cols] > 0.001).all(axis=1)].copy()
        log_message(log_file, f"After tip filter (>0.001): {len(filtered)} regions")

        if len(filtered) == 0:
            log_message(log_file, "No regions passed tip filter – empty outputs.")
            pd.DataFrame(columns=df.columns).to_csv(output_table, sep="\t", index=False)
            open(neutral_list, "w").close()
            return

        # ---- IQR filter on foreground relative branch length ----
        rb_col = f"{fg_branch}_rb"
        if rb_col not in filtered.columns:
            raise ValueError(f"Column {rb_col} missing. Columns: {list(filtered.columns)}")

        q1, q3 = filtered[rb_col].quantile([0.25, 0.75])
        neutral_mask = (filtered[rb_col] >= q1) & (filtered[rb_col] <= q3)
        neutralset = filtered[neutral_mask].copy()

        log_message(log_file, f"IQR filter ({rb_col}): Q1={q1:.6f}, Q3={q3:.6f} → {len(neutralset)} neutral regions")

        # ---- Write final outputs ----
        neutralset.to_csv(output_table, sep="\t", index=False)

        with open(neutral_list, "w") as f:
            f.write("\n".join(f"{r}.fa" for r in neutralset["region"]) + "\n")

        log_message(log_file, f"SUCCESS – neutral_table.txt ({len(neutralset)} rows), neutralset.txt")

    except Exception as e:
        log_message(log_file, f"FATAL ERROR: {e}\n{traceback.format_exc()}")
        sys.exit(1)


if __name__ == "__main__":
    main()

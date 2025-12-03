#!/usr/bin/env python3
"""
Select and filter pruned FASTA alignments.

Behavior:
- Recursively scan INPUT_DIR for FASTA files (*.fa, *.fasta, *.fa.prunned).
- For each file, require at least MIN_REAL_BASES A/C/G/T bases in every sequence.
- Keep up to MAX_KEEP files (deterministic first-N unless --random is set).
- Copy kept files into OUT_DIR (default: good_alignments/) preserving basename.
- Produce:
    * manifest file with absolute paths to kept files (--manifest)
    * goodalignments.txt listing kept files (one per line) (--good-list)
- Move failing files into TRASH_DIR (default: INPUT_DIR/TRASH).
- Exit code: 0 on success (even if no files kept), non-zero on fatal error.

Usage:
  scripts/select_and_filter.py --input pruned_neutral --min-bases 200 --max-keep 1000 \
      --out good_alignments --manifest selected_alignments/selected_alignments.list \
      --good-list goodalignments.txt
"""
import argparse
import os
import shutil
import random
from pathlib import Path
from typing import List
from Bio import SeqIO

VALID_BASES = set("ACGTacgt")

def count_real_bases(s: str) -> int:
    return sum(1 for c in s if c in VALID_BASES)

def file_is_valid(path: Path, min_real: int) -> bool:
    try:
        records = list(SeqIO.parse(str(path), "fasta"))
        if len(records) < 2:
            return False
        for rec in records:
            seq = str(rec.seq)
            # count real bases across whole sequence (not just leading)
            if count_real_bases(seq) < min_real:
                return False
        return True
    except Exception as e:
        # parsing error -> treat as invalid
        return False

def find_fasta_files(root: Path) -> List[Path]:
    patterns = ("**/*.fa", "**/*.fasta", "**/*.fa.prunned", "**/*.fa.good.fa")
    out = []
    for pat in patterns:
        out.extend([p for p in root.glob(pat) if p.is_file()])
    # de-duplicate and sort
    return sorted({p.resolve() for p in out})

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="pruned input directory")
    p.add_argument("--min-bases", type=int, default=200)
    p.add_argument("--max-keep", type=int, default=1000)
    p.add_argument("--out", default="good_alignments")
    p.add_argument("--manifest", default="selected_alignments/selected_alignments.list")
    p.add_argument("--good-list", default="goodalignments.txt")
    p.add_argument("--trash", default=None)
    p.add_argument("--random", action="store_true", help="sample randomly instead of first-N")
    args = p.parse_args()

    in_dir = Path(args.input)
    if not in_dir.exists():
        print(f"ERROR: input directory {in_dir} does not exist", flush=True)
        return 2

    trash_dir = Path(args.trash) if args.trash else in_dir / "TRASH"
    out_dir = Path(args.out)
    manifest = Path(args.manifest)
    good_list = Path(args.good_list)

    out_dir.mkdir(parents=True, exist_ok=True)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    trash_dir.mkdir(parents=True, exist_ok=True)

    candidates = find_fasta_files(in_dir)
    if not candidates:
        # write empty outputs and exit 0
        manifest.write_text("")
        good_list.write_text("")
        print("No candidate FASTA files found; exiting OK")
        return 0

    valid = []
    invalid = []

    print(f"Found {len(candidates)} candidate FASTA files under {in_dir}")
    for fp in candidates:
        # skip files already in output dir to avoid re-copying
        if fp.resolve().parent == out_dir.resolve():
            continue
        ok = file_is_valid(fp, args.min_bases)
        if ok:
            valid.append(fp)
        else:
            invalid.append(fp)

    if args.random:
        random.shuffle(valid)

    keep = valid[: args.max_keep]

    # copy kept files into out_dir, write manifest and good_list
    kept_paths = []
    for src in keep:
        dst = out_dir / src.name
        try:
            if not dst.exists():
                shutil.copy2(str(src), str(dst))
            kept_paths.append(str(dst.resolve()))
        except Exception as e:
            print(f"ERROR copying {src} -> {dst}: {e}")

    with manifest.open("w") as mf:
        for p in kept_paths:
            mf.write(p + "\n")

    with good_list.open("w") as gf:
        for p in kept_paths:
            gf.write(p + "\n")

    # move invalid files to trash (do not trash files that were kept or already in out_dir)
    for bad in invalid:
        try:
            if bad.resolve().parent == out_dir.resolve():
                continue
            dst = trash_dir / bad.name
            shutil.move(str(bad), str(dst))
        except Exception as e:
            print(f"WARNING: failed to move {bad} to trash: {e}")

    print(f"Kept {len(kept_paths)} files (written to {out_dir}); manifest: {manifest}")
    print(f"Moved {len(invalid)} invalid files to {trash_dir}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())

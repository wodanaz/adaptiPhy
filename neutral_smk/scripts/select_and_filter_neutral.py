#!/usr/bin/env python3
"""
Safe filter: copies good alignments from deep subdirectories.
Works correctly with msa_split --out-root pruned_neutral/{chrom}/{chrom}
"""
import argparse
import shutil
from pathlib import Path
from Bio import SeqIO

VALID_BASES = set("ACGTacgt")

def count_real_bases(seq: str) -> int:
    return sum(1 for c in seq if c in VALID_BASES)

def is_valid_fasta(path: Path, min_bases: int) -> bool:
    try:
        records = list(SeqIO.parse(str(path), "fasta"))
        if len(records) < 2:
            return False
        return all(count_real_bases(str(rec.seq)) >= min_bases for rec in records)
    except Exception:
        return False

def main():
    parser = argparse.ArgumentParser(description="Copy good alignments")
    parser.add_argument("--input", required=True, help="Top-level pruned directory (e.g. pruned_neutral/chr1)")
    parser.add_argument("--min-bases", type=int, default=200)
    parser.add_argument("--out", default="good_alignments")
    args = parser.parse_args()

    in_dir = Path(args.input)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not in_dir.exists():
        print(f"[ERROR] Input directory not found: {in_dir}")
        return 1

    candidates = []
    for pattern in ["**/*.fa", "**/*.fasta", "**/*.fa.prunned", "**/*.fa.good.fa"]:
        candidates.extend(in_dir.glob(pattern))

    candidates = sorted({p.resolve() for p in candidates if p.is_file()})
    print(f"[INFO] Found {len(candidates)} candidate FASTA files in {in_dir}")

    kept = 0
    for src in candidates:
        if src.resolve().parent == out_dir.resolve():
            continue
        if is_valid_fasta(src, args.min_bases):
            dst = out_dir / src.name
            if not dst.exists():
                shutil.copy2(src, dst)
            kept += 1

    print(f"[SUCCESS] Copied {kept} good alignments → {out_dir}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())

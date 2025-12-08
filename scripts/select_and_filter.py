#!/usr/bin/env python3
"""
Filter: copies good alignments from deep subdirectories to output directory.
Works correctly with msa_split --out-root query/{chrom}/{chrom}
"""
import argparse
import shutil
from pathlib import Path
from Bio import SeqIO

VALID_BASES = set("ACGTacgt")

def count_real_bases(seq: str) -> int:
    return sum(1 for c in seq if c in VALID_BASES)

def is_valid_fasta(path: Path, min_frac: float) -> bool:
    try:
        records = list(SeqIO.parse(str(path), "fasta"))
        if len(records) < 2:
            print(f"  [INVALID] {path.name}: Only {len(records)} sequences (<2 required)")
            return False
        # Check ALL sequences (not just ref) for sufficient ACGT fraction
        for i, rec in enumerate(records):
            seq = str(rec.seq)
            real_count = count_real_bases(seq)
            frac = real_count / len(seq) if len(seq) > 0 else 0
            if frac < min_frac:
                print(f"  [INVALID] {path.name}: Seq {i+1} ('{rec.id}') ACGT frac={frac:.3f} < {min_frac}")
                return False
        print(f"  [VALID] {path.name}: {len(records)} seqs, all fracs >= {min_frac}")
        return True
    except Exception as e:
        print(f"  [INVALID] {path.name}: Parse error - {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Copy good alignments")
    parser.add_argument("--input", required=True, help="Top-level pruned directory (e.g. query/chr1)")
    parser.add_argument("--min-frac", type=float, default=0.9)
    parser.add_argument("--out", default="good_alignments")
    args = parser.parse_args()
    in_dir = Path(args.input)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    if not in_dir.exists():
        print(f"[ERROR] Input directory not found: {in_dir}")
        return 1
    candidates = []
    patterns = ["**/*.fa", "**/*.fasta","**/*/*.fasta", "**/*/*.fa"] 
    for pattern in patterns:
        candidates.extend(in_dir.glob(pattern))
    candidates = sorted({p.resolve() for p in candidates if p.is_file()})
    print(f"[INFO] Found {len(candidates)} candidate FASTA files in {in_dir}")
    if not candidates:
        print(f"[WARNING] No FASTA files matched patterns: {patterns}")
        return 0
    copied = 0
    for src in candidates:
        if src.resolve().parent == out_dir.resolve():
            print(f"  [SKIP] {src.name}: Already in output")
            continue
        if is_valid_fasta(src, args.min_frac):
            dst = out_dir / src.name
            if not dst.exists():
                shutil.copy2(src, dst)
                print(f"  [COPIED] {src.name} -> {dst}")
            else:
                print(f"  [EXISTS] {src.name} already copied")
            copied += 1
    print(f"[SUCCESS] Copied {copied} good alignments → {out_dir}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())

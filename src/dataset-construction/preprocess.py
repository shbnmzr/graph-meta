#!/usr/bin/env python3
"""
Genome Cleaning & Deduplication for refs4 Dataset
---------------------------------------------------------------

This script processes your refs4 dataset structure:

    base_dir/
        train/prokaryote/
        train/plasmid/
        train/eukaryote/
        train/virus/
        val/...
        test/...

It writes cleaned genomes to a separate directory:

    clean_base_dir/
        train/prokaryote/*.cleaned.fna
        val/prokaryote/*.cleaned.fna
        ...

And performs:
    - Cleaning of each genome
    - Exact deduplication (MD5)
    - Mash-based near-identical deduplication
    - Optional selection of categories (--category)

"""

import argparse
import logging
import subprocess
import hashlib
from pathlib import Path
from Bio import SeqIO
import csv

# ----------------------------------------------------
# CONFIG
# ----------------------------------------------------

CATEGORIES = ["prokaryote", "plasmid", "eukaryote", "virus"]
AMBIGUOUS = set("NRYSWKMBVDH")

# ----------------------------------------------------
# FASTA UTILITIES
# ----------------------------------------------------

def safe_parse(path):
    """Safely parse FASTA records."""
    try:
        return list(SeqIO.parse(path, "fasta"))
    except Exception as e:
        logging.error(f"Failed to parse {path}: {e}")
        return []

def is_ambiguous(seq):
    """Count ambiguous bases in a sequence."""
    return sum(base in AMBIGUOUS for base in seq)

def compute_genome_hash(path):
    """Compute MD5 hash for entire genome (sorted concatenation of contigs)."""
    seqs = safe_parse(path)
    concat = "".join(sorted(str(rec.seq).upper() for rec in seqs))
    return hashlib.md5(concat.encode()).hexdigest()

# ----------------------------------------------------
# CLEANING
# ----------------------------------------------------

def clean_fasta(path, out_dir, min_len, max_ambig):
    """Clean a genome FASTA by removing short/ambiguous/duplicate contigs."""
    seen = set()
    kept = []
    before = 0

    for rec in safe_parse(path):
        before += 1
        seq = str(rec.seq).upper()

        if not seq:
            continue

        # Deduplicate contigs inside genome
        h = hashlib.md5(seq.encode()).hexdigest()
        if h in seen:
            continue
        seen.add(h)

        # Minimum contig length
        if len(seq) < min_len:
            continue

        # Ambiguity check
        if (is_ambiguous(seq) / len(seq)) > max_ambig:
            continue

        rec.seq = rec.seq.upper()
        kept.append(rec)

    if not kept:
        return None, before, 0

    out_path = out_dir / (path.stem + ".cleaned.fna")
    SeqIO.write(kept, out_path, "fasta")
    return out_path, before, len(kept)

# ----------------------------------------------------
# EXACT DEDUPLICATION (MD5)
# ----------------------------------------------------

def dedup_md5(clean_base, group, log_dir):
    logging.info(f"[EXACT DEDUP] Group={group}")

    splits = ["train", "val", "test"]
    genome_hashes = {s: {} for s in splits}

    # Compute hashes
    for split in splits:
        dir_path = clean_base / split / group
        if not dir_path.exists():
            continue
        for f in dir_path.glob("*.cleaned.fna"):
            genome_hashes[split][f] = compute_genome_hash(f)

    removed = []

    train_hashes = set(genome_hashes["train"].values())
    val_hashes   = set(genome_hashes["val"].values())

    # Remove val/test duplicates of train
    for split in ["val", "test"]:
        for f, h in list(genome_hashes[split].items()):
            if h in train_hashes and f.exists():
                logging.warning(f"Removing EXACT duplicate ({split}): {f.name}")
                f.unlink()
                removed.append((split, f.name, "EXACT_DUPLICATE_OF_TRAIN"))

    # Remove test duplicates of val
    for f, h in list(genome_hashes["test"].items()):
        if h in val_hashes and f.exists():
            logging.warning(f"Removing EXACT duplicate (test): {f.name}")
            f.unlink()
            removed.append(("test", f.name, "EXACT_DUPLICATE_OF_VAL"))

    # Report
    report = log_dir / f"{group}_exact_dedup.csv"
    with open(report, "w") as fh:
        w = csv.writer(fh)
        w.writerow(["split", "filename", "reason"])
        for row in removed:
            w.writerow(row)


# ----------------------------------------------------
# MASH NEAR-IDENTICAL DEDUPLICATION
# ----------------------------------------------------

def mash_sketch(paths, prefix):
    cmd = ["mash", "sketch", "-o", prefix] + [str(p) for p in paths]
    subprocess.run(cmd, check=True)
    return f"{prefix}.msh"

def mash_dist(query_msh, ref_msh):
    cmd = ["mash", "dist", ref_msh, query_msh]
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)

    distances = {}
    for line in res.stdout.strip().split("\n"):
        if not line.strip():
            continue
        ref_f, q_f, dist, _, _ = line.split("\t")
        distances[(Path(q_f).name, Path(ref_f).name)] = float(dist)
    return distances

def dedup_mash(clean_base, group, log_dir, threshold=0.05):
    logging.info(f"[MASH] Group={group}")

    def compare(ref_split, query_split):
        ref_dir = clean_base / ref_split / group
        q_dir   = clean_base / query_split / group

        ref_paths = list(ref_dir.glob("*.cleaned.fna"))
        q_paths   = list(q_dir.glob("*.cleaned.fna"))

        if not ref_paths or not q_paths:
            return []

        removed_local = []

        ref_sketch = mash_sketch(ref_paths, str(log_dir / f"{group}_{ref_split}"))

        for q in q_paths:
            q_sketch = mash_sketch([q], str(log_dir / f"{group}_{query_split}_{q.stem}"))
            distances = mash_dist(q_sketch, ref_sketch)

            for (qname, rname), d in distances.items():
                if d <= threshold:
                    logging.warning(
                        f"Removing NEAR IDENTICAL genome ({query_split}): "
                        f"{qname} ~ {rname} (d={d:.3f})"
                    )
                    q.unlink()
                    removed_local.append((query_split, qname, f"MASH_{d:.3f}"))
                    break

        return removed_local

    removed = []
    removed += compare("train", "val")
    removed += compare("train", "test")
    removed += compare("val", "test")

    report = log_dir / f"{group}_mash_dedup.csv"
    with open(report, "w") as fh:
        w = csv.writer(fh)
        w.writerow(["split", "filename", "reason"])
        for r in removed:
            w.writerow(r)


# ----------------------------------------------------
# MAIN
# ----------------------------------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--base_dir", required=True, help="Path to dataset/refs4/")
    p.add_argument("--clean_base_dir", required=True, help="Path to write cleaned genomes")
    p.add_argument("--log_dir", required=True)
    p.add_argument("--min_len", type=int, default=1000)
    p.add_argument("--max_ambig", type=float, default=0.05)
    p.add_argument("--min_contigs", type=int, default=1)

    # Run on selected categories only
    p.add_argument(
        "--category",
        nargs="*",
        choices=CATEGORIES,
        help="Optional: process only these categories. "
             "If omitted, all 4 categories are processed."
    )

    args = p.parse_args()

    base_dir = Path(args.base_dir)
    clean_base = Path(args.clean_base_dir)
    log_dir = Path(args.log_dir)

    log_dir.mkdir(parents=True, exist_ok=True)
    clean_base.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

    # Choose categories dynamically
    category_to_process = args.group if args.group else CATEGORIES
    logging.info(f"Processing category: {category_to_process}")

    # -----------------------------
    # MAIN LOOP
    # -----------------------------
    for group in category_to_process:
        logging.info(f"=== GROUP: {group.upper()} ===")

        # CLEANING
        for split in ["train", "val", "test"]:
            in_dir  = base_dir / split / group
            out_dir = clean_base / split / group
            out_dir.mkdir(parents=True, exist_ok=True)

            fasta_files = []
            for ext in ["*.fna", "*.fa", "*.fasta"]:
                fasta_files.extend(in_dir.glob(ext))

            for fasta in fasta_files:
                logging.info(f"[{split}/{group}] Cleaning: {fasta.name}")
                cleaned_path, before, after = clean_fasta(
                    fasta, out_dir, args.min_len, args.max_ambig
                )
                if cleaned_path and after < args.min_contigs and cleaned_path.exists():
                    cleaned_path.unlink()

        # EXACT DEDUP
        dedup_md5(clean_base, group, log_dir)

        # NEAR IDENTICAL DEDUP
        dedup_mash(clean_base, group, log_dir)

    logging.info("=== ALL DONE ===")


if __name__ == "__main__":
    main()

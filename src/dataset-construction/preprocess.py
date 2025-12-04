#!/usr/bin/env python3
"""
Genome Cleaning & Deduplication for refs4 Dataset
-------------------------------------------------

This script cleans & deduplicates each category, separately.

Valid categories (must match relabel step):
    prokaryote
    plasmid
    eukaryote
    viral

Input directory structure:
    base_dir/
        train/<category>/
        val/<category>/
        test/<category>/

Output directory structure:
    clean_base_dir/
        train/<category>/*.cleaned.fna
        val/<category>/*.cleaned.fna
        test/<category>/*.cleaned.fna

Steps performed:
  Cleaning: remove ambiguous, short, duplicate contigs
  Exact deduplication (MD5) across splits
  Mash-based near-identical removal across splits

Usage:

python3 preprocess.py \
    --base_dir dataset/refs4"\
    --clean_base_dir dataset/refs4_cleaned \
    --log_dir logs/dedup \
    --category prokaryote
"""

import argparse
import logging
import subprocess
import hashlib
from pathlib import Path
from Bio import SeqIO
import csv

# Allowed categories (MUST match relabel script output)
VALID_CATEGORIES = ["prokaryote", "plasmid", "eukaryote", "viral"]

AMBIGUOUS = set("NRYSWKMBVDH")

# ---------------------------
# FASTA Parsing Utils
# ---------------------------

def safe_parse(path):
    try:
        return list(SeqIO.parse(path, "fasta"))
    except Exception as e:
        logging.error(f"Failed to parse {path}: {e}")
        return []

def is_ambiguous(seq):
    return sum(b in AMBIGUOUS for b in seq)

# ---------------------------
# Cleaning
# ---------------------------

def clean_fasta(path, out_dir, min_len, max_ambig):
    """Clean FASTA: remove short, ambiguous, duplicate contigs."""
    seen = set()
    kept = []
    before = 0

    for rec in safe_parse(path):
        before += 1
        seq = str(rec.seq).upper()
        if not seq:
            continue

        # Deduplicate contigs
        h = hashlib.md5(seq.encode()).hexdigest()
        if h in seen:
            continue
        seen.add(h)

        # Filters
        if len(seq) < min_len:
            continue
        if (is_ambiguous(seq) / len(seq)) > max_ambig:
            continue

        rec.seq = rec.seq.upper()
        kept.append(rec)

    if not kept:
        return None, before, 0

    out_path = out_dir / (path.stem + ".cleaned.fna")
    SeqIO.write(kept, out_path, "fasta")
    return out_path, before, len(kept)

# ---------------------------
# Exact Deduplication (MD5)
# ---------------------------

def compute_genome_hash(path):
    seqs = safe_parse(path)
    concat = "".join(sorted(str(rec.seq).upper()) for rec in seqs)
    return hashlib.md5(concat.encode()).hexdigest()

def dedup_md5(clean_base, category, log_dir):
    logging.info(f"[MD5] Deduplicating category: {category}")

    splits = ["train", "val", "test"]
    hashes = {s: {} for s in splits}

    for split in splits:
        d = clean_base / split / category
        for f in d.glob("*.cleaned.fna"):
            hashes[split][f] = compute_genome_hash(f)

    removed = []
    train_hashes = set(hashes["train"].values())
    val_hashes = set(hashes["val"].values())

    # Remove val/test duplicates of train
    for split in ["val", "test"]:
        for f, h in list(hashes[split].items()):
            if h in train_hashes and f.exists():
                logging.warning(f"MD5 duplicate removed ({split}): {f.name}")
                f.unlink()
                removed.append([split, f.name, "MD5_DUPLICATE_OF_TRAIN"])

    # Remove test duplicates of val
    for f, h in list(hashes["test"].items()):
        if h in val_hashes and f.exists():
            logging.warning(f"MD5 duplicate removed (test): {f.name}")
            f.unlink()
            removed.append(["test", f.name, "MD5_DUPLICATE_OF_VAL"])

    report = log_dir / f"{category}_exact_dedup.csv"
    with open(report, "w") as fh:
        w = csv.writer(fh)
        w.writerow(["split", "filename", "reason"])
        for r in removed:
            w.writerow(r)

# ---------------------------
# Mash Near-Identical Dedupe
# ---------------------------

def mash_sketch(paths, prefix):
    cmd = ["mash", "sketch", "-o", prefix] + [str(p) for p in paths]
    subprocess.run(cmd, check=True)
    return f"{prefix}.msh"

def mash_dist(q_msh, r_msh):
    cmd = ["mash", "dist", r_msh, q_msh]
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    distances = {}
    for line in res.stdout.strip().split("\n"):
        if not line.strip():
            continue
        ref_f, q_f, d, _, _ = line.split("\t")
        distances[(Path(q_f).name, Path(ref_f).name)] = float(d)
    return distances

def dedup_mash(clean_base, category, log_dir, threshold=0.05):
    logging.info(f"[MASH] Deduplicating category: {category}")

    def compare(ref_split, query_split):
        ref_dir = clean_base / ref_split / category
        q_dir = clean_base / query_split / category

        ref_files = list(ref_dir.glob("*.cleaned.fna"))
        q_files = list(q_dir.glob("*.cleaned.fna"))

        if not ref_files or not q_files:
            return []

        removed_local = []
        ref_msh = mash_sketch(ref_files, str(log_dir / f"{category}_{ref_split}"))

        for q in q_files:
            q_msh = mash_sketch([q], str(log_dir / f"{category}_{query_split}_{q.stem}"))
            dists = mash_dist(q_msh, ref_msh)

            for (qname, rname), dist in dists.items():
                if dist <= threshold:
                    logging.warning(
                        f"Removing near-identical ({query_split}): "
                        f"{qname} ~ {rname} (d={dist:.3f})"
                    )
                    q.unlink()
                    removed_local.append([query_split, qname, f"MASH_{dist:.3f}"])
                    break
        return removed_local

    removed = []
    removed += compare("train", "val")
    removed += compare("train", "test")
    removed += compare("val", "test")

    report = log_dir / f"{category}_mash_dedup.csv"
    with open(report, "w") as fh:
        w = csv.writer(fh)
        w.writerow(["split", "filename", "reason"])
        for r in removed:
            w.writerow(r)

# ---------------------------
# MAIN
# ---------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base_dir", required=True, help="refs4 directory")
    parser.add_argument("--clean_base_dir", required=True, help="cleaned output directory")
    parser.add_argument("--log_dir", required=True)
    parser.add_argument("--category", required=True,
                        choices=VALID_CATEGORIES,
                        help="One of: prokaryote, plasmid, eukaryote, viral")
    parser.add_argument("--min_len", type=int, default=1000)
    parser.add_argument("--max_ambig", type=float, default=0.05)
    args = parser.parse_args()

    base_dir = Path(args.base_dir)
    clean_dir = Path(args.clean_base_dir)
    log_dir = Path(args.log_dir)
    category = args.category

    log_dir.mkdir(parents=True, exist_ok=True)
    clean_dir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

    logging.info(f"Processing category: {category}")

    # -------------------------
    # 1. CLEANING
    # -------------------------
    for split in ["train", "val", "test"]:
        in_dir = base_dir / split / category
        out_dir = clean_dir / split / category
        out_dir.mkdir(parents=True, exist_ok=True)

        fasta_files = []
        for ext in ["*.fna", "*.fa", "*.fasta"]:
            fasta_files.extend(in_dir.glob(ext))

        if not fasta_files:
            logging.warning(f"No FASTA files in {in_dir}")
            continue

        for fasta in fasta_files:
            logging.info(f"[{split}/{category}] Cleaning {fasta.name}")
            cleaned_path, before, after = clean_fasta(
                fasta, out_dir, args.min_len, args.max_ambig
            )
            if cleaned_path and after == 0 and cleaned_path.exists():
                cleaned_path.unlink()

    # -------------------------
    # 2. EXACT DEDUP (MD5)
    # -------------------------
    dedup_md5(clean_dir, category, log_dir)

    # -------------------------
    # 3. NEAR IDENTICAL (MASH)
    # -------------------------
    dedup_mash(clean_dir, category, log_dir)

    logging.info("Finished processing category: " + category)


if __name__ == "__main__":
    main()

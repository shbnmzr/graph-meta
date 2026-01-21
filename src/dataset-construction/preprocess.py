#!/usr/bin/env python3
"""
Perform splitting, classification, cleaning, and deduplication.

Workflow:
  1. Iterate splits in hierarchical order (Train -> Val -> Test).
  2. Classify records (Prokaryote/Eukaryote/Viral/Plasmid).
  3. Filter: Remove sequences with high ambiguity or short length.
  4. Clean: Strip non-ACGT characters.
  5. Dedup (Exact): Check MD5 against ALL previously seen sequences (preventing data leakage).
  6. Write: Save to final directory and update manifest.
  7. Post-process: Run Mash (MinHash) to remove near-identical sequences across splits.

Usage:
  python preprocess.py \
      --base_dir ./raw_data \
      --out_dir ./processed_data \
      --category bacteria \
      --mash_threshold 0.05
"""

import argparse
import csv
import gzip
import hashlib
import json
import logging
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Union, Optional

from Bio import SeqIO
from Bio.Seq import Seq

# ----------------------------------
# CONFIGURATION
# ----------------------------------

CATEGORY_TO_DOMAIN = {
    "archaea": "prokaryote",
    "bacteria": "prokaryote",
    "fungi": "eukaryote",
    "protozoa": "eukaryote",
    "virus": "viral",
    "viral": "viral",
}

# Regex pre-compilation
REGEX_NON_ACGT = re.compile(r"[^ACGT]")
REGEX_SANITIZER = re.compile(r"[^A-Za-z0-9._-]+")


# ----------------------------------
# UTILS
# ----------------------------------

def sanitize_id(s: str) -> str:
    """Sanitize string for safe filenames."""
    return REGEX_SANITIZER.sub("_", s or "")


def get_replicon_type(description: str) -> str:
    """Classify replicon based on description keywords."""
    d = (description or "").lower()
    if any(w in d for w in ["phage", "virus", "viral", "bacteriophage"]):
        return "viral"
    if "plasmid" in d:
        return "plasmid"
    if any(w in d for w in ["chromosome", "complete genome", "chromosomal"]):
        return "chromosomal"
    return "unknown"


def get_class4(category: str, description: str) -> str:
    """Determine the high-level class (prokaryote, eukaryote, viral, plasmid)."""
    rtype = get_replicon_type(description)
    if rtype in {'plasmid', 'viral'}:
        return rtype
    return CATEGORY_TO_DOMAIN.get(category.lower(), 'unknown')


def iter_fasta(path: Path):
    """Robust iterator for plain or gzipped FASTA."""
    open_func = gzip.open if path.suffix == ".gz" else open
    try:
        with open_func(path, "rt", encoding="utf-8", errors="ignore") as fh:
            yield from SeqIO.parse(fh, "fasta")
    except Exception as e:
        logging.error(f"Error reading {path.name}: {e}")


# ----------------------------------
# CORE PROCESSOR
# ----------------------------------

class BioProcessor:
    def __init__(self, args):
        self.args = args
        self.seen_md5_global: Set[str] = set()

        # Stats tracking
        self.stats = defaultdict(int)

        # Host map storage
        self.host_map = defaultdict(lambda: {"chromosome_accessions": [], "plasmids": {}})

        # Setup directories
        self.meta_dir = self.args.out_dir / "metadata"
        self.meta_dir.mkdir(parents=True, exist_ok=True)

    def is_high_quality(self, seq_str: str) -> bool:
        """
        Check quality BEFORE cleaning.
        If a sequence is 50% Ns, we reject it rather than stripping it to half size.
        """
        if not seq_str:
            return False

        # 1. Ambiguity Check
        # Count non-ACGT characters (like N, R, Y)
        # Note: We assume input is uppercase for efficiency
        ambiguous_count = sum(1 for b in seq_str if b not in "ACGT")
        if (ambiguous_count / len(seq_str)) > self.args.max_ambig:
            return False

        return True

    def process_sequence(self, rec) -> Optional[Seq]:
        """
        Cleans and hashes sequence. Returns None if duplicate or low quality.
        """
        seq_upper = str(rec.seq).upper()

        # 1. Intra-file duplicate check (simple caching handled in loop) & Quality
        if not self.is_high_quality(seq_upper):
            self.stats["skipped_low_quality"] += 1
            return None

        # 2. Clean (Strip non-ACGT)
        cleaned_seq_str = REGEX_NON_ACGT.sub("", seq_upper)

        # 3. Length Check (Post-cleaning)
        if len(cleaned_seq_str) < self.args.min_len:
            self.stats["skipped_short"] += 1
            return None

        # 4. Global Deduplication (MD5)
        # We prevent data leakage by ensuring no sequence in Val/Test
        # is identical to one seen in Train.
        seq_hash = hashlib.md5(cleaned_seq_str.encode("utf-8")).hexdigest()

        if seq_hash in self.seen_md5_global:
            self.stats["skipped_duplicate_md5"] += 1
            return None

        self.seen_md5_global.add(seq_hash)

        return Seq(cleaned_seq_str)

    def extract_plasmid_name(self, desc: str, rec_id: str) -> str:
        """Helper to extract plasmid name from description."""
        d_lower = desc.lower()
        if "plasmid" in d_lower:
            parts = d_lower.split("plasmid", 1)
            if len(parts) > 1 and parts[1].strip():
                # Recover original casing using index
                idx = d_lower.find("plasmid") + 7
                return desc[idx:].strip()
        return rec_id

    def run_mash_dedup(self, category_name: str):
        """
        Runs Mash (MinHash) to remove near-identical sequences across splits.
        Requires 'mash' to be installed in the system path.
        """
        if shutil.which("mash") is None:
            logging.warning("MASH not found in PATH. Skipping near-identical deduplication.")
            return

        logging.info("Starting Mash deduplication...")

        # Define comparisons: Train vs Val, Train vs Test, Val vs Test
        pairs = [("train", "val"), ("train", "test"), ("val", "test")]

        for ref_split, query_split in pairs:
            ref_dir = self.args.out_dir / ref_split / category_name
            q_dir = self.args.out_dir / query_split / category_name

            ref_files = list(ref_dir.glob("*.fna"))
            q_files = list(q_dir.glob("*.fna"))

            if not ref_files or not q_files:
                continue

            # Sketch Reference
            sketch_prefix = str(self.meta_dir / f"mash_{ref_split}")
            subprocess.run(
                ["mash", "sketch", "-o", sketch_prefix] + [str(p) for p in ref_files],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )
            ref_msh = f"{sketch_prefix}.msh"

            # Check Queries
            for q_file in q_files:
                # Sketch single query
                q_sketch = str(self.meta_dir / f"temp_q")
                subprocess.run(
                    ["mash", "sketch", "-o", q_sketch, str(q_file)],
                    check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
                )

                # Dist
                res = subprocess.run(
                    ["mash", "dist", ref_msh, f"{q_sketch}.msh"],
                    capture_output=True, text=True
                )

                # Parse Results
                for line in res.stdout.strip().split("\n"):
                    if not line: continue
                    parts = line.split("\t")
                    dist = float(parts[2])

                    if dist <= self.args.mash_threshold:
                        logging.warning(f"Mash Dup ({dist:.4f}): Removing {q_file.name} (matches {parts[0]})")
                        q_file.unlink()
                        self.stats["skipped_mash_duplicate"] += 1
                        break  # File removed, stop checking against other refs

    def run(self):
        manifest_path = self.meta_dir / f"{self.args.category}_manifest.csv"

        with open(manifest_path, "w", newline="") as mf:
            writer = csv.DictWriter(mf, fieldnames=[
                "split", "category", "class4", "replicon_type",
                "host_assembly", "accession", "description", "path", "source_file"
            ])
            writer.writeheader()

            # STRICT ORDER: Train -> Val -> Test
            # This ensures Val/Test are deduped against Train, not vice versa.
            for split in ["train", "val", "test"]:
                input_split_dir = self.args.base_dir / split / self.args.category

                if not input_split_dir.exists():
                    logging.warning(f"Skipping missing directory: {input_split_dir}")
                    continue

                files = sorted(list(input_split_dir.glob("*")))
                logging.info(f"Processing {split} ({len(files)} files)...")

                for fpath in files:
                    if fpath.name.startswith("."): continue  # skip hidden

                    assembly_id = sanitize_id(fpath.stem)

                    # Track local hashes to dedup contigs within the same file
                    # (Though global hash set also covers this, this is explicit)

                    for rec in iter_fasta(fpath):
                        desc = rec.description or rec.id
                        rtype = get_replicon_type(desc)
                        cls = get_class4(self.args.category, desc)

                        if cls == "unknown" and not self.args.keep_unknown:
                            continue

                        # PROCESS SEQUENCE (Clean & Dedup)
                        processed_seq = self.process_sequence(rec)
                        if processed_seq is None:
                            continue

                        rec.seq = processed_seq

                        # HOST MAPPING
                        if rtype == "plasmid":
                            label = self.extract_plasmid_name(desc, rec.id)
                            self.host_map[assembly_id]["plasmids"][rec.id] = label
                        elif rtype == "chromosomal":
                            if rec.id not in self.host_map[assembly_id]["chromosome_accessions"]:
                                self.host_map[assembly_id]["chromosome_accessions"].append(rec.id)

                        # WRITE TO DISK
                        out_cls_dir = self.args.out_dir / split / cls
                        out_cls_dir.mkdir(parents=True, exist_ok=True)

                        rec_id = sanitize_id(rec.id)
                        out_name = f"{fpath.stem}__{rec_id}.fna"
                        out_path = out_cls_dir / out_name

                        SeqIO.write([rec], out_path, "fasta")
                        self.stats["records_written"] += 1

                        writer.writerow({
                            "split": split, "category": self.args.category,
                            "class4": cls, "replicon_type": rtype,
                            "host_assembly": assembly_id, "accession": rec.id,
                            "description": desc, "path": str(out_path),
                            "source_file": str(fpath)
                        })

            # Save Host Map
            with open(self.meta_dir / f"{self.args.category}_host_map.json", "w") as jf:
                json.dump(self.host_map, jf, indent=2)

            # Run Mash Deduplication (Post-processing)
            self.run_mash_dedup(self.args.category)

            logging.info("Processing Complete.")
            logging.info(f"Stats: {json.dumps(self.stats, indent=2)}")


# ----------------------------------
# MAIN ENTRY POINT
# ----------------------------------

def main():
    ap = argparse.ArgumentParser(description="Genome ETL: Split, Clean, Dedup (MD5+Mash)")
    ap.add_argument("--base_dir", required=True, type=Path, help="Input raw directory")
    ap.add_argument("--out_dir", required=True, type=Path, help="Output clean directory")
    ap.add_argument("--category", required=True, help="Taxonomic category (bacteria, viral, etc.)")

    # Filter Params
    ap.add_argument("--min_len", type=int, default=1000, help="Min length after cleaning")
    ap.add_argument("--max_ambig", type=float, default=0.05, help="Max fraction of Ns allowed (0.0-1.0)")
    ap.add_argument("--mash_threshold", type=float, default=0.05, help="Mash distance threshold for dedup")
    ap.add_argument("--keep_unknown", action="store_true", help="Keep records classified as unknown")

    args = ap.parse_args()

    # Logging Setup
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.StreamHandler()]
    )

    # Validate Category
    cat = args.category.lower()
    if cat == "virus": cat = "viral"
    if cat not in CATEGORY_TO_DOMAIN:
        logging.error(f"Invalid category: {cat}. Must be one of {list(CATEGORY_TO_DOMAIN.keys())}")
        return
    args.category = cat

    # Run Pipeline
    processor = BioProcessor(args)
    processor.run()


if __name__ == "__main__":
    main()

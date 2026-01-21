#!/usr/bin/env python3
"""
Split, clean & relabel to 4 classes, with host–plasmid mapping.

Classes: prokaryote, eukaryote, plasmid, viral

Pipeline:
  1. Reads <base_dir>/{train,val,test}/<category>/*.(fasta|gz)
  2. Classifies & Cleans sequences (keeps only ACGT)
  3. Writes per-record FASTAs to <out_dir>
  4. Generates a metadata manifest CSV
  5. Generates a Host-Plasmid JSON map
"""

import argparse
import csv
import gzip
import json
import logging
import re
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Optional, Union, Dict

from Bio import SeqIO
from Bio.Seq import Seq

# ----------------------------------
# CONFIG
# ----------------------------------

CATEGORY_TO_DOMAIN = {
    "archaea": "prokaryote",
    "bacteria": "prokaryote",
    "fungi": "eukaryote",
    "protozoa": "eukaryote",
    "virus": "viral",
    "viral": "viral",
}


# ----------------------------------
# HELPERS
# ----------------------------------

def clean_sequence(seq: Union[str, Seq]) -> str:
    """
    Uppercase and strip ambiguous bases. Keep only A,C,G,T.
    Returns a string for Regex performance.
    """
    return re.sub(r"[^ACGT]", "", str(seq).upper())


def replicon_type(description: str) -> str:
    """Classify replicon based on description keywords."""
    d = (description or "").lower()

    if any(w in d for w in ["phage", "virus", "viral", "bacteriophage"]):
        return "viral"

    if "plasmid" in d:
        return "plasmid"

    if any(w in d for w in ["chromosome", "complete genome", "chromosomal"]):
        return "chromosomal"

    return "unknown"


def class4_for_record(category: str, description: str) -> str:
    """Map replicon to high-level classes, respecting plasmid/viral overrides."""
    rtype = replicon_type(description)

    # 1. Description overrides category (e.g., plasmid in bacteria folder is plasmid)
    if rtype in {'plasmid', 'viral'}:
        return rtype

    # 2. Category default
    return CATEGORY_TO_DOMAIN.get(category.lower(), 'unknown')


def sanitize_id(s: str) -> str:
    """Sanitize string for use in filenames."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s or "")


def iter_fasta_records(path: Path):
    """Yields SeqRecords from plain or gzipped files."""
    open_func = gzip.open if str(path).endswith(".gz") else open
    with open_func(path, "rt", encoding="utf-8", errors="ignore") as fh:
        yield from SeqIO.parse(fh, "fasta")


# ----------------------------------
# PROCESSING LOGIC
# ----------------------------------

def process_file(
        fpath: Path,
        split: str,
        category: str,
        out_dir: Path,
        writer: csv.DictWriter,
        host_map: Dict,
        keep_unknown: bool,
) -> int:
    written = 0
    assembly_id = sanitize_id(fpath.stem)

    try:
        for rec in iter_fasta_records(fpath):
            desc = rec.description or rec.id
            rtype = replicon_type(desc)
            cls = class4_for_record(category, desc)

            # Filter unknowns
            if cls == "unknown" and not keep_unknown:
                continue

            # Clean Sequence (Str -> Regex -> Check -> Seq)
            clean_seq_str = clean_sequence(rec.seq)
            if len(clean_seq_str) < 1000:
                continue

            # Update record with clean sequence
            rec.seq = Seq(clean_seq_str)

            # Host-plasmid mapping
            if rtype == "plasmid":
                plasmid_label = rec.id
                d_lower = desc.lower()
                if "plasmid" in d_lower:
                    # Attempt to extract plasmid name
                    parts = d_lower.split("plasmid", 1)
                    if len(parts) > 1 and parts[1].strip():
                        # Use original case by index logic
                        idx = d_lower.find("plasmid") + 7
                        plasmid_label = desc[idx:].strip()

                host_map[assembly_id]["plasmids"][rec.id] = plasmid_label

            elif rtype == "chromosomal":
                # Only add if not already present
                if rec.id not in host_map[assembly_id]["chromosome_accessions"]:
                    host_map[assembly_id]["chromosome_accessions"].append(rec.id)

            # Determine Output Paths
            out_cls_dir = out_dir / split / cls
            rec_id = sanitize_id(rec.id)
            out_name = f"{fpath.stem}__{rec_id}.fna"
            out_path = out_cls_dir / out_name

            # Write to disk
            out_cls_dir.mkdir(parents=True, exist_ok=True)
            SeqIO.write([rec], out_path, "fasta")

            # Write Manifest
            writer.writerow({
                "split": split,
                "category": category,
                "class4": cls,
                "replicon_type": rtype,
                "host_assembly": assembly_id,
                "accession": rec.id,
                "description": desc,
                "path": str(out_path),
                "source_file": str(fpath),
            })

            written += 1

    except Exception as e:
        logging.exception(f"Failed parsing {fpath.name}: {e}")

    return written


def process_category(
        base_dir: Path,
        out_dir: Path,
        category: str,
        writer: csv.DictWriter,
        host_map: Dict,
        keep_unknown: bool,
        delete_source: bool,
        trash_dir: Optional[Path],
):
    """Iterates over train/val/test splits for a given category."""

    for split in ["train", "val", "test"]:
        in_dir = base_dir / split / category

        if not in_dir.exists():
            logging.warning(f"Skipping missing dir: {in_dir}")
            continue

        files = sorted([
            p for p in in_dir.glob("*")
            if p.suffix in {'.fa', '.fna', '.fasta', '.gz'}
        ])

        logging.info(f"[{split}/{category}] Processing {len(files)} files...")

        for fpath in files:
            n = process_file(
                fpath, split, category, out_dir, writer, host_map, keep_unknown
            )

            # Cleanup Source Files
            if n > 0:
                if trash_dir:
                    dst = trash_dir / fpath.name
                    trash_dir.mkdir(parents=True, exist_ok=True)
                    shutil.move(str(fpath), str(dst))
                elif delete_source:
                    fpath.unlink()
            elif n == 0:
                logging.debug(f"Skipped cleanup for {fpath.name} (no records written)")


# ----------------------------------
# MAIN
# ----------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base_dir", required=True, type=Path)
    ap.add_argument("--out_dir", required=True, type=Path)
    ap.add_argument("--category", required=True)
    ap.add_argument("--keep_unknown", action="store_true")
    ap.add_argument("--delete_source", action="store_true")
    ap.add_argument("--trash_dir", type=Path, default=None)
    ap.add_argument("--log_level", default="INFO")
    args = ap.parse_args()

    logging.basicConfig(level=args.log_level, format="%(asctime)s [%(levelname)s] %(message)s")

    # Validate Category
    cat = args.category.lower()
    if cat == "virus": cat = "viral"
    if cat not in CATEGORY_TO_DOMAIN:
        logging.error(f"Category '{cat}' not recognized in configuration.")
        return

    # Initialize Global Host Map
    host_map = defaultdict(lambda: {"chromosome_accessions": [], "plasmids": {}})

    # Setup Manifest
    meta_root = args.out_dir.parent if args.out_dir.name == "refs4" else args.out_dir
    manifest_path = meta_root / "metadata" / f"{cat}_refs4_manifest.csv"
    host_map_path = meta_root / "metadata" / f"{cat}_host_plasmids.json"

    # Open Manifest Writer
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    with manifest_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            "split", "category", "class4", "replicon_type",
            "host_assembly", "accession", "description", "path", "source_file"
        ])
        writer.writeheader()

        process_category(
            args.base_dir, args.out_dir, cat, writer, host_map,
            args.keep_unknown, args.delete_source, args.trash_dir
        )

    # Write Host Map
    host_map_path.parent.mkdir(parents=True, exist_ok=True)
    with host_map_path.open("w") as f:
        json.dump(host_map, f, indent=2)


if __name__ == "__main__":
    main()

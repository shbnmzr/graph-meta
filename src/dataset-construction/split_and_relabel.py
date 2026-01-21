#!/usr/bin/env python3
"""
Split, clean & relabel to 4 classes (with host–plasmid mapping)

Classes: prokaryote, eukaryote, plasmid, viral

Reads:
  <base_dir>/{train,val,test}/<category>/*.(fna|fa|fasta|*.gz)

Writes per-record FASTAs:
  <out_dir>/{train,val,test}/{prokaryote|eukaryote|plasmid|viral}/

Writes manifest:
  <meta_root>/metadata/<category>_refs4_manifest.csv
  where meta_root = out_dir/.. if out_dir name == "refs4", else out_dir

Also writes host–plasmid map:
  <meta_root>/metadata/<category>_host_plasmids.json
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
from typing import List, Optional, Union
from Bio import SeqIO
from Bio.Seq import Seq

# ----------------------------------
# CONFIG / CONSTANTS
# ----------------------------------

# Map taxonomic category → high-level class for chromosomal replicons
CATEGORY_TO_DOMAIN = {
    "archaea": "prokaryote",
    "bacteria": "prokaryote",
    "fungi": "eukaryote",
    "protozoa": "eukaryote",
    "virus": "viral",
    "viral": "viral",
}

# Store host↔plasmid/co-chromosome relationships
HOST_MAP = defaultdict(lambda: {
    "chromosome_accessions": [],
    "plasmids": {}
})


# ----------------------------------
# HELPERS
# ----------------------------------

def clean_sequence(seq: Union[str, Seq]) -> Union[str, Seq]:
    """Uppercase and strip ambiguous bases. Keep only A,C,G,T."""
    seq_str = str(seq)
    cleaned = re.sub(r'[^ACGT]', '', seq_str.upper())

    if isinstance(seq, Seq):
        return Seq(cleaned)

    return cleaned

def replicon_type(description: str) -> str:
    """
    Classify a replicon within an assembly:
      - 'plasmid'
      - 'viral' (phage, virus, viral, bacteriophage)
      - 'chromosomal'
      - 'unknown'
    """
    d = (description or "").lower()

    # viral replicons
    if any(w in d for w in ["phage", "virus", "viral", "bacteriophage"]):
        return "viral"

    # plasmids
    if "plasmid" in d:
        return "plasmid"

    # chromosomes / complete genomes
    if any(w in d for w in ["chromosome", "complete genome", "chromosomal"]):
        return "chromosomal"

    return "unknown"


def class4_for_record(category: str, description: str) -> str:
    """
    Map a replicon to the 4 high-level classes:
      - prokaryote
      - eukaryote
      - plasmid
      - viral
      - unknown
    """
    rtype = replicon_type(description)

    # Check description-based overrides
    if rtype in {'plasmid', 'viral'}:
        return rtype

    # Check category-based mapping
    # return 'unknown'  if missing
    return CATEGORY_TO_DOMAIN.get(category.lower(), 'unknown')


def sanitize_id(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s or "")


def iter_fasta_records(path: Path):
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as fh:
            yield from SeqIO.parse(fh, "fasta")
    else:
        with open(path, "rt", encoding="utf-8", errors="ignore") as fh:
            yield from SeqIO.parse(fh, "fasta")


def find_fasta_files(root: Path) -> List[Path]:
    exts = ["*.fna", "*.fa", "*.fasta", "*.fna.gz", "*.fa.gz", "*.fasta.gz"]
    files: List[Path] = []
    for ext in exts:
        files.extend(root.glob(ext))
    return sorted(files)


def open_manifest(manifest_path: Path):
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    fh = manifest_path.open("w", newline="")
    writer = csv.DictWriter(
        fh,
        fieldnames=[
            "split",
            "category",
            "class4",
            "replicon_type",
            "host_assembly",
            "accession",
            "description",
            "path",
            "source_file",
        ],
    )
    writer.writeheader()
    return fh, writer


# ----------------------------------
# CORE PROCESSING
# ----------------------------------

def process_file(
    fpath: Path,
    split: str,
    category: str,
    out_dir: Path,
    writer: csv.DictWriter,
    keep_unknown: bool,
    dry_run: bool,
) -> int:

    written = 0
    assembly_id = sanitize_id(fpath.stem)  # e.g. GCA_000005825.2

    try:
        for rec in iter_fasta_records(fpath):
            desc = rec.description or rec.id
            rtype = replicon_type(desc)
            cls = class4_for_record(category, desc)

            if cls == "unknown" and not keep_unknown:
                continue

            # Cleaned sequence
            clean_seq = clean_sequence(rec.seq)
            if len(clean_seq) < 1000:
                continue

            # Host-plasmid mapping
            if rtype == "plasmid":
                d_lower = desc.lower()
                if "plasmid" in d_lower:
                    label_part = desc.split("plasmid", 1)[-1].strip()
                    plasmid_label = label_part if label_part else rec.id
                else:
                    plasmid_label = rec.id

                HOST_MAP[assembly_id]["plasmids"][rec.id] = plasmid_label

            elif rtype == "chromosomal":
                if rec.id not in HOST_MAP[assembly_id]["chromosome_accessions"]:
                    HOST_MAP[assembly_id]["chromosome_accessions"].append(rec.id)

            # Output directory
            out_cls_dir = out_dir / split / cls
            if not dry_run:
                out_cls_dir.mkdir(parents=True, exist_ok=True)

            rec_id = sanitize_id(rec.id)
            out_name = f"{fpath.stem}__{rec_id}.fna"
            out_path = out_cls_dir / out_name

            if not dry_run:
                SeqIO.write([rec], str(out_path), "fasta")

            # Write manifest entry
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
        logging.exception(f"Failed parsing {fpath}: {e}")

    return written


def delete_or_trash(fpath: Path, trash_dir: Optional[Path], dry_run: bool):
    if trash_dir is not None:
        trash_dir.mkdir(parents=True, exist_ok=True)
        dst = trash_dir / fpath.name
        if dry_run:
            logging.info(f"[DRY] Would move source → {dst}")
        else:
            shutil.move(str(fpath), str(dst))
    else:
        if not dry_run:
            fpath.unlink(missing_ok=True)


def process_split_category(
    base_dir: Path,
    out_dir: Path,
    split: str,
    category: str,
    writer: csv.DictWriter,
    keep_unknown: bool,
    delete_source: bool,
    trash_dir: Optional[Path],
    dry_run: bool,
):
    in_dir = base_dir / split / category
    if not in_dir.exists():
        logging.warning(f"[{split}/{category}] directory not found: {in_dir}")
        return

    fasta_files = find_fasta_files(in_dir)
    logging.info(f"[{split}/{category}] Found {len(fasta_files)} files")

    for fpath in fasta_files:
        logging.info(f"Processing {fpath}")
        n = process_file(fpath, split, category, out_dir, writer, keep_unknown, dry_run)

        if n > 0 and (delete_source or trash_dir is not None):
            delete_or_trash(fpath, trash_dir, dry_run)
        elif n == 0:
            logging.warning(f"No records written from {fpath}; skipping deletion.")


# ----------------------------------
# CLI
# ----------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Split, clean, and relabel FASTA records into 4 classes; optionally delete originals and record host–plasmid relationships."
    )
    ap.add_argument("--base_dir", required=True, type=Path)
    ap.add_argument("--out_dir",  required=True, type=Path)
    ap.add_argument("--category", required=True,
                    help="archaea | bacteria | fungi | protozoa | virus | viral")
    ap.add_argument("--keep_unknown", action="store_true")
    ap.add_argument("--delete_source", action="store_true")
    ap.add_argument("--trash_dir", type=Path, default=None)
    ap.add_argument("--dry_run", action="store_true")
    ap.add_argument("--log_level", default="INFO",
                    choices=["DEBUG","INFO","WARNING","ERROR"])
    args = ap.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s"
    )

    cat = args.category.lower()
    if cat not in {"archaea", "bacteria", "fungi", "protozoa", "virus", "viral"}:
        raise SystemExit(f"Unsupported category: {args.category}")

    # Normalize category name
    if cat == "virus":
        cat = "viral"

    # Manifest path
    meta_root = args.out_dir.parent if args.out_dir.name == "refs4" else args.out_dir

    manifest_path = meta_root / "metadata" / f"{cat}_refs4_manifest.csv"
    host_map_path = meta_root / "metadata" / f"{cat}_host_plasmids.json"

    # Open manifest unless dry run
    if args.dry_run:
        fh = None
        writer = None
    else:
        fh, writer = open_manifest(manifest_path)

    try:
        for split in ["train", "val", "test"]:
            if args.dry_run:
                class _Dummy:
                    def writerow(self, row): pass
                writer = _Dummy()

            process_split_category(
                base_dir=args.base_dir,
                out_dir=args.out_dir,
                split=split,
                category=cat,
                writer=writer,
                keep_unknown=args.keep_unknown,
                delete_source=args.delete_source,
                trash_dir=args.trash_dir,
                dry_run=args.dry_run,
            )
    finally:
        if not args.dry_run and fh is not None:
            fh.close()

        # Write host–plasmid map
        if not args.dry_run:
            host_map_path.parent.mkdir(parents=True, exist_ok=True)
            with host_map_path.open("w") as f:
                json.dump(HOST_MAP, f, indent=2)


if __name__ == "__main__":
    main()

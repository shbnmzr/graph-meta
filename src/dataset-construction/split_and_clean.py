#!/usr/bin/env python3
"""
Split, clean & relabel to 4 classes, with plasmid-host mapping

Classes: prokaryote, eukaryote, plasmid, virus

Reads:
  <base_dir>/{train,val,test}/<category>/*.(fna|fa|fasta|*.gz)

Writes per-record FASTAs:
  <out_dir>/{train,val,test}/{prokaryote|eukaryote|plasmid|virus}/

Writes manifest:
  <meta_root>/metadata/<category>_refs4_manifest.csv
  where meta_root = out_dir/.. if out_dir name == "refs4", else out_dir

Writes host–plasmid map:
  <meta_root>/metadata/<category>_host_plasmids.json

Flags:
  --delete_source    delete original file after successful split (>=1 record written)
  --trash_dir DIR    move original files here instead of deleting (overrides --delete_source)
  --dry_run          do not write or delete/move anything; just log actions

Usage:
  python split_and_relabel.py \
    --base_dir /path/to/my_dataset \
    --out_dir  /path/to/my_dataset/refs4 \
    --category bacteria \
    --delete_source
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
from typing import List, Optional
from Bio import SeqIO
from Bio.Seq import Seq

# ----------------------------------
# CONFIG / CONSTANTS
# ----------------------------------

CATEGORY_TO_DOMAIN = {
    "archaea": "prokaryote",
    "bacteria": "prokaryote",
    "fungi": "eukaryote",
    "protozoa": "eukaryote",
    "virus": "virus",
}

# host assembly → { "chromosome_accessions": [...], "plasmids": {plasmid_acc: description} }
HOST_MAP = defaultdict(lambda: {
    "chromosome_accessions": [],
    "plasmids": {}
})


# ----------------------------------
# HELPERS
# ----------------------------------

def clean_sequence(seq: str) -> str:
    """Uppercase and strip ambiguous bases. Keep only A,C,G,T."""
    seq = seq.upper()
    return "".join(b for b in seq if b in "ACGT")


def replicon_type(description: str) -> str:
    """
    Classify a replicon within an assembly:
      - 'plasmid'
      - 'virus' (phage, viral, bacteriophage)
      - 'chromosomal'
      - 'unknown'
    """
    d = (description or "").lower()

    # phage / viral
    if any(w in d for w in ["phage", "virus", "viral", "bacteriophage"]):
        return "virus"

    # plasmids
    if "plasmid" in d:
        return "plasmid"

    # chromosomes / complete genomes
    if any(w in d for w in ["chromosome", "complete genome", "chromosomal"]):
        return "chromosomal"

    # default
    return "unknown"


def class4_for_record(category: str, description: str) -> str:
    """
    Map a replicon to the 4 high-level classes:
      - prokaryote
      - eukaryote
      - plasmid
      - virus
      - unknown
    """
    rtype = replicon_type(description)
    if rtype == "plasmid":
        return "plasmid"
    if rtype == "virus":
        return "virus"

    dom = CATEGORY_TO_DOMAIN.get(category.lower())
    if dom in {"prokaryote", "eukaryote", "virus"}:
        return dom
    return "unknown"


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
    """
    Process a single source FASTA file.
    Returns number of records successfully written from this source file.
    """
    written = 0
    assembly_id = sanitize_id(fpath.stem)  # e.g. GCA_000005825.2

    try:
        for rec in iter_fasta_records(fpath):
            desc = rec.description or rec.id
            rtype = replicon_type(desc)
            cls = class4_for_record(category, desc)

            if cls == "unknown" and not keep_unknown:
                continue

            # Clean sequence: keep only A,C,G,T and drop short ones
            clean_seq = clean_sequence(str(rec.seq))
            if len(clean_seq) < 1000:
                continue

            # Update SeqRecord with cleaned sequence
            rec.seq = Seq(clean_seq)

            # Track host–plasmid relationships
            if rtype == "plasmid":
                # Try to extract plasmid label after the word "plasmid"
                d_lower = desc.lower()
                if "plasmid" in d_lower:
                    # crude split; safe enough for later lookup
                    label_part = desc.split("plasmid", 1)[-1].strip()
                    plasmid_label = label_part if label_part else rec.id
                else:
                    plasmid_label = rec.id

                HOST_MAP[assembly_id]["plasmids"][rec.id] = plasmid_label

            elif rtype == "chromosomal":
                if rec.id not in HOST_MAP[assembly_id]["chromosome_accessions"]:
                    HOST_MAP[assembly_id]["chromosome_accessions"].append(rec.id)

            # Prepare output directory
            out_cls_dir = out_dir / split / cls
            if not dry_run:
                out_cls_dir.mkdir(parents=True, exist_ok=True)

            src_stem = fpath.stem
            rec_id = sanitize_id(rec.id)
            out_name = f"{src_stem}__{rec_id}.fna"
            out_path = out_cls_dir / out_name

            if dry_run:
                logging.info(f"[DRY] Would write {out_path}")
            else:
                SeqIO.write([rec], str(out_path), "fasta")

            # Write manifest row
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
            logging.info(f"Moved source → {dst}")
    else:
        if dry_run:
            logging.info(f"[DRY] Would DELETE source {fpath}")
        else:
            fpath.unlink(missing_ok=True)
            logging.info(f"Deleted source {fpath}")


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
) -> None:
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
    ap.add_argument("--base_dir", required=True, type=Path,
                    help="Base dir with train/val/test/<category>/*")
    ap.add_argument("--out_dir",  required=True, type=Path,
                    help="Output dir for refs4/<split>/<class4>/")
    ap.add_argument("--category", required=True,
                    help="archaea | bacteria | fungi | protozoa | virus")
    ap.add_argument("--keep_unknown", action="store_true",
                    help="Keep records that cannot be classified")
    ap.add_argument("--delete_source", action="store_true",
                    help="Delete original file after successful split")
    ap.add_argument("--trash_dir", type=Path, default=None,
                    help="Move originals here instead of deleting")
    ap.add_argument("--dry_run", action="store_true",
                    help="Log actions without writing or deleting")
    ap.add_argument("--log_level", default="INFO",
                    choices=["DEBUG","INFO","WARNING","ERROR"])
    args = ap.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s"
    )

    cat = args.category.lower()
    if cat not in {"archaea","bacteria","fungi","protozoa","virus"}:
        raise SystemExit(f"Unsupported category: {args.category}")

    # Decide metadata root
    meta_root = args.out_dir.parent if args.out_dir.name == "refs4" else args.out_dir

    manifest_path = meta_root / "metadata" / f"{cat}_refs4_manifest.csv"
    host_map_path = meta_root / "metadata" / f"{cat}_host_plasmids.json"

    if args.dry_run:
        logging.info(f"[DRY] Would write manifest to {manifest_path}")
        logging.info(f"[DRY] Would write host–plasmid map to {host_map_path}")

    # Open manifest (unless dry run)
    if args.dry_run:
        fh = None
        writer = None
    else:
        fh, writer = open_manifest(manifest_path)

    try:
        for split in ["train","val","test"]:
            if args.dry_run:
                # create a temporary no-op writer that just logs
                class _Dummy:
                    def writerow(self, row): logging.debug(f"[DRY] manifest row: {row}")
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
            logging.info(f"Wrote manifest: {manifest_path}")

        # Write host–plasmid map
        if not args.dry_run:
            host_map_path.parent.mkdir(parents=True, exist_ok=True)
            with host_map_path.open("w") as f:
                json.dump(HOST_MAP, f, indent=2)
            logging.info(f"Wrote host–plasmid map: {host_map_path}")


if __name__ == "__main__":
    main()

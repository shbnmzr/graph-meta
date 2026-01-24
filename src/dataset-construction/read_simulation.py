import os
import shutil
from pathlib import Path
import logging
from tqdm.notebook import tqdm
from typing import List
import math
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np


def combine_references(split_dir: Path, output_fasta: Path) -> bool:
    """
    Combines all the .fna files in a split into one fasta file
    """
    logging.info(f'[1/4] Combining references from {split_dir}...')

    with open(output_fasta, 'w') as outfile:
        files = list(split_dir.glob('**/*.fna'))
        if not files:
            logging.warning(f'No files found in {split_dir}')
            return False

        logging.info(f'Found {len(files)} genomes.')
        with open(output_fasta, 'wb') as outfile:
            for fname in tqdm(files, desc='Merging Genomes', unit='file'):
                with open(fname, 'rb') as infile:
                    shutil.copyfileobj(infile, outfile)

                    outfile.write(b'\n')


        return True


def run_simulation(split: str,
                   input_dir: Path,
                   output_dir: Path,
                   N_READS: int,
                   BATCH_SIZE: int,
                   DISTRIBUTION_MODE: str):
    logging.info(f'--- Processing Split: {split.upper()} ---')
    logging.info(f'[Scientific Mode: {DISTRIBUTION_MODE}] Simulating {N_READS} reads...')

    all_genomes = list(input_dir.rglob('*.fna'))
    if not all_genomes: all_genomes = list(input_dir.rglob('*.fasta'))

    if not all_genomes:
        logging.warning(f'No genomes found in {input_dir}')
        return

    # Batching Setup
    num_batches = math.ceil(len(all_genomes) / BATCH_SIZE)
    reads_per_batch = max(1, int(N_READS / num_batches))

    final_r1 = output_dir / 'simulated_R1.fastq'
    final_r2 = output_dir / 'simulated_R2.fastq'

    if final_r1.exists(): final_r1.unlink()
    if final_r2.exists(): final_r2.unlink()

    temp_dir = Path("/content/temp_scientific_work")
    if temp_dir.exists(): shutil.rmtree(temp_dir)
    temp_dir.mkdir(parents=True, exist_ok=True)

    batch_iterator = get_file_batches(all_genomes, BATCH_SIZE)

    for i, batch_files in enumerate(tqdm(batch_iterator, total=num_batches, desc=f'Simulating {split}')):

        # 1. Prepare Batch File
        batch_ref = temp_dir / 'batch_ref.fasta'
        merge_fasta_batch(batch_files, batch_ref)

        # 2. Parse Records into Memory (CRITICAL FIX)
        # We must convert the File Path -> List of SeqRecords here
        try:
            with open(batch_ref) as f:
                batch_records = list(SeqIO.parse(f, "fasta"))
        except:
            continue

        if not batch_records: continue

        # 3. Calculate Distribution
        read_counts = calculate_read_counts(batch_records, reads_per_batch, mode=DISTRIBUTION_MODE)

        # 4. Generate Reads
        batch_r1 = temp_dir / f'batch_{i}_R1.fastq'
        batch_r2 = temp_dir / f'batch_{i}_R2.fastq'

        with open(batch_r1, 'w') as f1, open(batch_r2, 'w') as f2:
            simulate_reads(batch_records, f1, f2, read_counts)

        # 5. Append & Cleanup
        if batch_r1.exists() and batch_r2.exists():
            append_fastq(batch_r1, final_r1)
            append_fastq(batch_r2, final_r2)

        for f in temp_dir.iterdir(): f.unlink()

    if temp_dir.exists(): temp_dir.rmdir()

    logging.info(f'Simulation completed on {split}')
    logging.info(f'Outputs: {final_r1}, {final_r2}')


def append_fastq(source: Path, dest: Path) -> None:
    """
    Appends source FASTQ content to destination FASTQ.
    """
    with open(dest, 'ab') as outfile:
        with open(source, 'rb') as infile:
            shutil.copyfileobj(infile, outfile)


def merge_fasta_batch(files: List, output_path: Path) -> None:
    """
    Merges a small list of files into one temporary fasta
    """
    with open(output_path, 'w') as outfile:
        for fname in files:
            file_stem = fname.stem
            with open(fname, 'r') as infile:
                for line in infile:
                    if line.startswith('>'):
                        outfile.write(f'>{file_stem}__{line.strip()[1:]}\n')
                    else:
                        outfile.write(line)

            outfile.write('/n')


def get_file_batches(files: List, batch_size: int):
    """
    A generator that yields successive chinks of files
    """
    for i in range(0, len(files), batch_size):
        yield files[i: i + batch_size]

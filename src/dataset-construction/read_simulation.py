import shutil
import logging
import math
import random
import numpy as np
from pathlib import Path
from typing import List, Generator, Dict
from Bio import SeqIO
from tqdm.notebook import tqdm

# Configure Logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')


def get_file_batches(files: List[Path],
                     batch_size: int) -> Generator[List[Path], None, None]:
    """
    Generator that yields successive chunks of files.
    """
    for i in range(0, len(files), batch_size):
        yield files[i: i + batch_size]


def merge_fasta_batch(files: List[Path], output_path: Path) -> None:
    """
    Merges a list of files into one temp fasta, renaming headers to ensure uniqueness.
    Optimized to minimize memory usage during merge.
    """
    with open(output_path, 'w') as outfile:
        for fname in files:
            file_stem = fname.stem
            try:
                with open(fname, 'r') as infile:
                    for line in infile:
                        if line.startswith('>'):
                            # Efficient string formatting
                            clean_header = line.strip()[1:]
                            outfile.write(f'>{file_stem}__{clean_header}\n')
                        else:
                            outfile.write(line)
                outfile.write('\n')
            except Exception as e:
                logging.warning(f"Skipping corrupt file {fname}: {e}")


def calculate_read_counts(records: List,
                          total_reads_target: int,
                          mode: str = 'uniform') -> Dict[str, int]:
    """
    Distributes reads among records using vectorized numpy operations for speed.
    """
    # 1. Vectorize Lengths
    lengths = np.array([len(rec.seq) for rec in records])

    if mode == 'uniform':
        weights = lengths.astype(float)
    elif mode == 'lognormal':
        # Log-Normal (Mean=0, Sigma=1.0) is standard for metagenomic abundance
        abundances = np.random.lognormal(mean=0, sigma=1.0, size=len(lengths))
        weights = lengths * abundances
    else:
        raise ValueError(f"Unknown mode: {mode}")

    # 2. Normalize & Distribute
    probs = weights / weights.sum()
    read_counts = np.random.multinomial(total_reads_target, probs)

    return {rec.id: count for rec, count in zip(records, read_counts)}


def simulate_reads(records: List,
                   r1_handle,
                   r2_handle,
                   read_counts_map: Dict[str, int],
                   insert_size: int,
                   read_len: int) -> None:
    """
    Shreds sequences into reads.
    Optimized: Calculates invariant strings outside the loop.
    """
    qual_string = "I" * read_len

    for record in records:
        target_reads = read_counts_map[record.id]
        seq = record.seq
        seq_len = len(seq)

        # Skip fragments shorter than insert size + buffer
        if seq_len < insert_size + 50:
            continue

        for _ in range(target_reads):
            # 1. Select coordinates
            start = random.randint(0, seq_len - insert_size)
            end = start + insert_size

            # 2. Extract fragment
            fragment = seq[start:end]

            # 3. Generate Pairs
            read1_seq = fragment[:read_len]
            read2_seq = fragment[-read_len:].reverse_complement()

            # 4. Write to disk
            header = f"@{record.id}_{start}"

            r1_handle.write(f"{header}/1\n{read1_seq}\n+\n{qual_string}\n")
            r2_handle.write(f"{header}/2\n{read2_seq}\n+\n{qual_string}\n")


def append_fastq(source: Path, dest: Path) -> None:
    """
    Appends source FASTQ content to destination FASTQ using binary stream (fastest copy).
    """
    if not source.exists(): return
    with open(dest, 'ab') as outfile:
        with open(source, 'rb') as infile:
            shutil.copyfileobj(infile, outfile)


def run_simulation(split: str,
                   input_dir: Path,
                   output_dir: Path,
                   n_reads: int,
                   batch_size: int,
                   distribution_mode: str,
                   read_len: int = 150,
                   insert_size: int = 400):
    logging.info(f'--- Processing Split: {split.upper()} ---')
    logging.info(f'[Mode: {distribution_mode}] Simulating {n_reads} reads...')

    # 1. Find Genomes
    all_genomes = list(input_dir.rglob('*.fna'))
    if not all_genomes:
        all_genomes = list(input_dir.rglob('*.fasta'))

    if not all_genomes:
        logging.warning(f'No genomes found in {input_dir}')
        return

    # 2. Setup Batching
    num_batches = math.ceil(len(all_genomes) / batch_size)
    reads_per_batch = max(1, int(n_reads / num_batches))

    # 3. Setup Outputs
    final_r1 = output_dir / 'simulated_R1.fastq'
    final_r2 = output_dir / 'simulated_R2.fastq'
    if final_r1.exists(): final_r1.unlink()
    if final_r2.exists(): final_r2.unlink()

    # 4. Temp Directory (Local SSD for speed)
    temp_dir = Path("/content/temp_scientific_work")
    if temp_dir.exists(): shutil.rmtree(temp_dir)
    temp_dir.mkdir(parents=True, exist_ok=True)

    # 5. Execution Loop
    batch_iterator = get_file_batches(all_genomes, batch_size)

    try:
        for i, batch_files in enumerate(tqdm(batch_iterator, total=num_batches, desc=f'Simulating {split}')):

            # A. Merge Batch
            batch_ref = temp_dir / 'batch_ref.fasta'
            merge_fasta_batch(batch_files, batch_ref)

            # B. Parse Records
            try:
                # list() loads batch into RAM
                with open(batch_ref) as f:
                    batch_records = list(SeqIO.parse(f, "fasta"))
            except Exception as e:
                logging.warning(f"Failed to parse batch {i}: {e}")
                continue

            if not batch_records: continue

            # C. Calculate Counts
            read_counts = calculate_read_counts(batch_records, reads_per_batch, mode=distribution_mode)

            # D. Generate Reads
            batch_r1 = temp_dir / f'batch_{i}_R1.fastq'
            batch_r2 = temp_dir / f'batch_{i}_R2.fastq'

            with open(batch_r1, 'w') as f1, open(batch_r2, 'w') as f2:
                simulate_reads(batch_records, f1, f2, read_counts, insert_size, read_len)

            # E. Append to Final
            append_fastq(batch_r1, final_r1)
            append_fastq(batch_r2, final_r2)

            # F. Cleanup Batch
            for f in temp_dir.iterdir():
                f.unlink()

    finally:
        # cleanup even if error occurs
        if temp_dir.exists(): shutil.rmtree(temp_dir)

    logging.info(f'Simulation completed on {split}')
    logging.info(f'Outputs: {final_r1}, {final_r2}')

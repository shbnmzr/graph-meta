import os
import shutil
from pathlib import Path
import logging
from tqdm.notebook import tqdm
from typing import Tuple


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


def read_simulation(combined_ref, out_dir: Path, N_READS: int, CORES: int) -> Tuple[Path, Path]:
    """
    Runs InSilicoSeq to generate paired-end reads
    """
    logging.info(f'[2 / 4] Simulating {N_READS} reads with InSilicoSeq...')

    out_prefix = out_dir / 'simulated'

    cmd = (f'iss generate --genome {combined_ref}'
           f'--model miseq --n_reads {N_READS}'
           f'--output {out_prefix} --cpus {CORES}')

    if os.system(cmd) != 0:
        raise RuntimeError('InSilicoSeq failed to generate reads')

    return out_prefix.with_name('simulated_R1.fastq'), out_prefix.with_name('simulated_R2.fastq')

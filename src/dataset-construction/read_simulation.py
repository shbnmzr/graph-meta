import os
import shutil
from pathlib import Path
import glob
import logging
from tqdm.notebook import tqdm


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

import os
import shutil
import logging
import pandas as pd
from pathlib import Path
from tqdm.notebook import tqdm

MEGAHIT_CMD_TEMPLATE = (
    'megahit -1 {r1} -2 {r2} -o {out_dir} '
    '--k-min 21 --k-max 121 --k-step 20 '
    '--min-contig-len 500 --out-prefix final'
)

MINIMAP_CMD_TEMPLATE = 'minimap2 -x asm5 {ref} {query} > {paf_out} 2>/dev/null'

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

def create_combined_reference(split_dir: Path, output_fasta: Path) -> bool:
    """
    Merges all the reference genomes in a split into a single fasta file for Minimap2
    :param split_dir: Path to the directory where genomes are located
    :param output_fasta: Path to the output fasta file
    :return:
    """
    logging.info(f'Creating temporary combined reference for labeling...')

    files = list(split_dir.glob('*.fna'))
    if not files: files = list(split_dir.glob('*.fasta'))

    if not files:
        logging.error('No reference genomes found!')
        return False

    with open(output_fasta, 'w') as outfile:
        for fname in tqdm(files, desc='Merging Refs'):
            file_stem = fname.stem
            try:
                with open(fname, 'r') as infile:
                    for line in infile:
                        if line.startswith('>'):
                            clean_header = line.strip()[1:]
                            outfile.write(f'>{file_stem}__{clean_header}\n')
                        else:
                            outfile.write(line)
                outfile.write('\n')
            except:
                continue
        return True

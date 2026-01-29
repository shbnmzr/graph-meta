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


def run_assembly(r1: Path, r2: Path, out_dir: Path):
    """
    Creates the De Novo Assembly graph using MEGAHIT
    :param r1:
    :param r2:
    :param out_dir:
    :return:
    """
    if out_dir.exists():
        logging.info('Assembly directory already exists, removing to start fresh...')
        shutil.rmtree(out_dir)

    logging.info('Running MEGAHIT...')
    cmd = MEGAHIT_CMD_TEMPLATE.format(r1=r1, r2=r2, out_dir=out_dir)

    exit_code = os.system(cmd)
    if exit_code != 0:
        raise RuntimeError(f'MEGAHIT failed with exit code {exit_code}')

    contigs = out_dir / 'final.contigs.fa'
    graph = out_dir / 'final.gfa'

    if not contigs.exists():
        raise FileNotFoundError('MEGAHIT finished but no contigs files found')

    return contigs, graph


def run_labeling(contigs: Path, ref_fasta: Path, output_csv: Path):
    """
    Maps contigs to reference genomes to establish the ground truth
    """
    paf_temp = output_csv.with_suffix('.paf')

    # Run minimap2
    logging.info('Mapping contigs to references using Minimap2...')
    cmd = MINIMAP_CMD_TEMPLATE.format(ref=ref_fasta, query=contigs, paf_out=paf_temp)
    os.system(cmd)

    # Parse PAF to CSV
    logging.info('Parsing mapping results...')

    labels = []
    if paf_temp.exists():
        with open(paf_temp, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 6: continue

                contig_id = parts[0]
                ref_header = parts[5]

                labels.append({
                    'contig_id': contig_id,
                    'ref_header': ref_header,
                    'length': int(parts[1])
                })
        df = pd.DataFrame(labels)

        df = df.sort_values(by='length', ascending=False).drop_duplicates('contig_id')
        df.to_csv(output_csv, index=False)

        logging.info(f'Labels saved to {output_csv} ({len(df)} contigs labeled)')

        paf_temp.unlink()

    else:
        logging.warning('Minimap2 failed to generate output.')

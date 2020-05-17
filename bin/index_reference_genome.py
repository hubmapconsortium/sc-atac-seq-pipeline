#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from shutil import copy
from subprocess import run
from typing import List, Optional, Tuple

from utils import find_base_index_path, normalize_whitespace

SNAPTOOLS_INDEX_REFERENCE_COMMAND = [
    'snaptools',
    'index-genome',
    '--aligner=bwa',
    # from Dockerfile
    '--path-to-aligner=/usr/local/bin',
    '--input-fasta={reference_genome}',
]
SAMTOOLS_FAIDX_COMMAND = [
    'samtools',
    'faidx',
    '{reference_genome}',
]

EXTRACT_ALIGNMENT_INDEX_COMMAND = [
    'tar',
    '-xf',
    '{alignment_index}',
]

def index_reference(reference_genome: Path):
    index_dest_dir = Path('index')
    index_dest_dir.mkdir()
    snaptools_command = [
        piece.format(reference_genome=reference_genome)
        for piece in SNAPTOOLS_INDEX_REFERENCE_COMMAND
    ]
    print('Running', ' '.join(snaptools_command))
    run(snaptools_command, check=True, cwd=index_dest_dir)

    samtools_command = [
        piece.format(reference_genome=reference_genome)
        for piece in SAMTOOLS_FAIDX_COMMAND
    ]
    print('Running', ' '.join(snaptools_command))
    run(samtools_command)

def extract_alignment_index(alignment_index: Path) -> Path:
    index_dest_dir = Path('index')
    index_dest_dir.mkdir()
    command = [
        piece.format(alignment_index=alignment_index)
        for piece in EXTRACT_ALIGNMENT_INDEX_COMMAND
    ]
    print('Running', ' '.join(command))
    run(command, check=True, cwd=index_dest_dir)
    return index_dest_dir

def read_alignment_index_seq_lengths(bwa_ann_file: Path) -> List[Tuple[str, int]]:
    # Returning a list instead of a generator so we can be nice
    # and show the length in the calling function
    seq_lengths = []
    with open(bwa_ann_file) as f:
        for line in f:
            pieces = line.strip().split()
            if pieces[0] != '0' or len(pieces) <= 3:
                continue
            sequence_name = pieces[1]
            if pieces[4] == 'bp':
                length = int(pieces[3])
            elif pieces[4].startswith('LN:'):
                length = int(pieces[4].split(':')[1])
            seq_lengths.append((sequence_name, length))
    return seq_lengths

def create_size_index_from_alignment_index(alignment_index_dir: Path) -> Path:
    """
    Parses supplementary data in the BWA index to create a sequence length
    index in the same format as `samtools faidx`. This is *NOT* the same
    as what `samtools` would produce, and is not usable for everything that
    a "real" .fai file would be.

    Specifically, we only know the lengths of each sub-sequence, not:
     * offset in bytes in the reference FASTA file
     * number of base pairs in each line of a sequence
     * number of bytes in each line of a sequence
    so we write 0 values for these. This is is specifically useful for
    `snaptools snap-pre`, so this temporary file should not be exposed as
    an output of the overall workflow.

    :param alignment_index_dir: Directory containing a BWA index
    :return: Path to fai-like file
    """
    base_index_path = find_base_index_path(alignment_index_dir)
    bwa_ann_file = base_index_path.parent / f'{base_index_path.name}.ann'
    seq_lengths = read_alignment_index_seq_lengths(bwa_ann_file)
    print('Read lengths for', len(seq_lengths), 'sequences from', bwa_ann_file)

    # no concrete benefit to using Path for a filename in the current directory,
    # but I like being consistent and using Paths for everything filesystem-related
    output_file = Path(f'{base_index_path.stem}.fai')
    with open(output_file, 'w') as f:
        # This just needs to work in a quick-and-dirty way, so not using
        # anything more structured like pandas or even the Python CSV module
        for sequence_name, length in seq_lengths:
            line_data = [sequence_name, str(length)] + ['0'] * 3
            print('\t'.join(line_data), file=f)

    return output_file

def index_if_necessary(
        reference_genome: Optional[Path],
        alignment_index: Optional[Path],
        size_index: Optional[Path],
):
    if reference_genome:
        if alignment_index or size_index:
            message = """
            Found reference genome and [alignment index or size index]. Providing a
            reference genome in FASTA format is only necessary if it should be
            indexed by this pipeline; either omit the reference genome or omit
            the precomputed index files.
            """
            raise ValueError(normalize_whitespace(message))
        print('Indexing reference genome')
        index_reference(reference_genome)

    if size_index and not alignment_index:
        raise ValueError('Alignment index is required if providing size index.')

    if alignment_index:
        index_dir = extract_alignment_index(alignment_index)
        if size_index:
            size_index_dest = Path() / size_index.name
            print('Copying', size_index, 'to', size_index_dest)
            copy(size_index, size_index_dest)
        else:
            message = """
            Found alignment index, but no size index (from `samtools faidx`)
            provided. Creating temporary size index from alignment index.
            """
            print(normalize_whitespace(message))
            size_index = create_size_index_from_alignment_index(index_dir)
            print('Wrote dummy genome size index to:', size_index)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('--reference-genome', type=Path)

    indexed_ref_group = p.add_argument_group(
        title='pre-indexed reference genome',
        description='Reference genome index arguments',
    )
    indexed_ref_group.add_argument(
        '--alignment-index',
        type=Path,
        help=normalize_whitespace(
            """
            Index to reference genome, stored as a tar archive. Will be
            extracted with `tar -xf`, so the archive can be compressed
            with any method understood by `tar` (e.g. '.gz', '.bz2',
            '.xz').
            """
        )
    )
    indexed_ref_group.add_argument(
        '--size-index',
        type=Path,
        help='Size index of reference genome, as created by `samtools faidx`.',
    )
    args = p.parse_args()

    index_if_necessary(args.reference_genome, args.alignment_index, args.size_index)

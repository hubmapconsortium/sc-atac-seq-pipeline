#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from shutil import copy
from subprocess import run
from typing import Optional

from utils import normalize_whitespace

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

def extract_alignment_index(alignment_index: Path):
    index_dest_dir = Path('index')
    index_dest_dir.mkdir()
    command = [
        piece.format(alignment_index=alignment_index)
        for piece in EXTRACT_ALIGNMENT_INDEX_COMMAND
    ]
    print('Running', ' '.join(command))
    run(command, check=True, cwd=index_dest_dir)

def index_if_necessary(
        reference_genome: Optional[Path],
        alignment_index: Optional[Path],
        size_index: Optional[Path],
):
    if reference_genome:
        if alignment_index or size_index:
            message = """
            Found reference genome and alignment index or size index. Providing a
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
        extract_alignment_index(alignment_index)
        if not size_index:
            message = """
            Size index (from `samtools faidx`) is required if providing alignment
            index. Conversion of BWA index supplementary data to the output format
            of `samtools faidx` is not yet implemented. 
            """
            raise NotImplementedError(normalize_whitespace(message))
        size_index_dest = Path() / size_index.name
        print('Copying', size_index, 'to', size_index_dest)
        copy(size_index, size_index_dest)

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

#!/usr/bin/env python3
from argparse import ArgumentParser
from os import fspath
from pathlib import Path
from subprocess import run
from typing import Iterable, Optional

DEFAULT_ALIGNMENT_INDEX_PATH = Path('/opt/supplementary-data/bwa-index')

SNAPTOOLS_ALIGN_COMMAND = [
    'snaptools',
    'align-paired-end',
    '--input-reference',
    '{base_index_path}',
]

def find_base_index_path(alignment_index: Path) -> str:
    """
    :param alignment_index: Directory to BWA alignment index
    :return: "Base" path to alignment reference. This is a `str` and
    not a `Path` because it doesn't actually exist on disk
    """
    # Could use any of the files expected to be part of the BWA index,
    # and could check for the presence of all of them, but may as well
    # just use one
    bwt_files = list(alignment_index.glob('*.bwt'))
    assert len(bwt_files) == 1
    return fspath(alignment_index / bwt_files[0].stem)

def align_reads(
        alignment_index: Optional[Path],
        extra_snaptools_args: Iterable[str],
):
    if alignment_index is None:
        base_index_path = find_base_index_path(DEFAULT_ALIGNMENT_INDEX_PATH)
    else:
        base_index_path = find_base_index_path(alignment_index)

    print('Alignment index base path:', base_index_path)
    print('Extra snaptools args:', extra_snaptools_args)

    snaptools_command = [
        piece.format(base_index_path=base_index_path)
        for piece in SNAPTOOLS_ALIGN_COMMAND
    ]
    snaptools_command.extend(extra_snaptools_args)
    print('Running', ' '.join(snaptools_command))
    run(snaptools_command, check=True)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('--alignment-index', type=Path, nargs='?')
    p.add_argument(
        'snaptools_align_arg',
        nargs='*',
        help='Extra arguments passed directly to `snaptools align-paired-end`.'
    )
    args = p.parse_args()

    align_reads(args.alignment_index, args.snaptools_align_arg)

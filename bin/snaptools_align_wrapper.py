#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from typing import Iterable, Optional

DEFAULT_ALIGNMENT_INDEX_PATH = Path('/opt/supplementary-data/bwa-index')

def align_reads(
        alignment_index: Optional[Path],
        extra_snaptools_args: Iterable[str],
):
    print('Alignment index:', alignment_index)
    print('Extra snaptools args:', extra_snaptools_args)

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
